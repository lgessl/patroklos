#' @title Build a fitter with integrated cross-validation from a fitter
#' @description This function operator takes a patroklos-compliant fitter and 
#' builds a patroklos-compliant fitter with integrated cross validation from 
#' it.
#' @param fitter A patroklos-compliant fitter.
#' @param error A string or a function. How to measure the prediction 
#' error. If a character string, it must be one of "error_rate", "neg_roc_auc", or 
#' "neg_binomial_log_likelihood". If a function, it must take two arguments: `y` 
#' and `y_hat`, and return a numeric scalar. The returned fitter will use it to 
#' calculate the goodness of validated predictions. 
#' @param select logical. If `TRUE`, fitter will only return the model that 
#' minimizes the error. If `FALSE`, fitter will return a `ptk_hypertune` object
#' with all models and their errors.
#' @return A patroklos-compliant fitter with integrated cross-validation tuning.
#' @export
hypertune <- function(
    fitter, 
    error = c("error_rate", "neg_roc_auc", "neg_binomial_log_likelihood"),
    select = FALSE
){
    if (is.character(error)) {
        error <- match.arg(error)
        error_fun <- paste0("get_", error)
    } else if (is.function(error)) {
        error_fun <- error
    } else {
        stop("`error` must be a character or a function.")
    }
    stopifnot(is.logical(select))
    force(fitter)
    force(error_fun)

    cv_fitter <- function(x, y_bin, y_cox, ...) {
        # Get all possible combinations of hyperparameters
        grid <- expand.grid(list(...), stringsAsFactors = FALSE)
        fit_obj_list <- lapply(seq(nrow(grid)), function(i) {
            do.call(fitter, c(list(x = x, y_bin = y_bin, y_cox = y_cox), 
                as.list(grid[i, ])))
        })
        # Get rid of NAs, i.e. invalid hyperparameters
        na_bool <- is.na(fit_obj_list)
        fit_obj_list <- fit_obj_list[!na_bool]
        grid <- grid[!na_bool, ]
        # Build ptk_hypertune S3 object
        ptk_hypertune <- list()
        # Promote validated predictions to top level
        ptk_hypertune$val_predict_list <- lapply(fit_obj_list, function(x) x$val_predict)
        ptk_hypertune$lambda <- apply(grid, 1, function(r) paste(names(grid), r, 
            sep = "=", collapse = ", "))
        # Evaluate according to error
        error_v <- vapply(ptk_hypertune$val_predict_list, function(y_hat) {
            y_yhat <- intersect_by_names(y_bin, y_hat, rm_na = c(TRUE, FALSE))
            do.call(error_fun, y_yhat)
        }, numeric(1))
        ptk_hypertune$lambda_min_index <- which.min(error_v)
        ptk_hypertune$lambda_min <- ptk_hypertune$lambda[ptk_hypertune$lambda_min_index]
        ptk_hypertune$val_predict <- ptk_hypertune$val_predict_list[[ptk_hypertune$lambda_min_index]]
        ptk_hypertune$val_error <- error_v
        ptk_hypertune$min_error <- error_v[ptk_hypertune$lambda_min_index]
        ptk_hypertune$error_name <- error_fun
        class(ptk_hypertune) <- "ptk_hypertune"
        if (select) {
            fit_obj_list[-ptk_hypertune$lambda_min_index] <- NA
            ptk_hypertune$val_predict_list <- NA
        }
        ptk_hypertune$fit_obj_list <- fit_obj_list
        ptk_hypertune
    }
}

#' @export
predict.ptk_hypertune <- function(
    object,
    newx,
    lambda_index = object$lambda_min_index,
    ...
){
    obj <- object$fit_obj_list[[lambda_index]]
    if (length(obj) == 1 && is.na(obj)) 
        stop("Desired model is NA. Did accidentally you set `select = TRUE`?")
    predict(object = obj, newx = newx, ...)
}

# A decorator to tune those hyperparameters one can tune for any fitter (or 
# for a *uni*versal fitter). Right now, this concerns 
# combine_n_max_categorical_features; time_cutoffs should follow.
unitune <- function(fitter) {

    unifitter <- function(x, y_bin, y_cox, combine_n_max_categorical_features, ...) {
        core <- function(n_max_combo) {
            x <- trim_combos(x, n_max_combo) 
            fit_obj <- fitter(x, y_bin, y_cox, ...)
            fit_obj$combine_n_max_categorical_features <- n_max_combo
            fit_obj
        }
        fit_obj_list <- lapply(combine_n_max_categorical_features, core)
        fit_obj_list[[which.min(sapply(fit_obj_list, function(x)
            x$min_error))]]
    }
}

# Remove columns with combinations of comprising more than
# combine_n_max_categorical_features
trim_combos <- function(x, combine_n_max_categorical_features) {
    keep <- sapply(stringr::str_split(colnames(x), "&"), length) <= 
        combine_n_max_categorical_features
    x_slim <- x[, keep]
    attr(x_slim, "li_var_suffix") <- attr(x, "li_var_suffix")
    x_slim
}