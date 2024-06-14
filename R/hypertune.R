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
#' @return A patroklos-compliant fitter with integrated cross-validation tuning.
#' @export
hypertune <- function(
    fitter, 
    error = c("error_rate", "neg_roc_auc", "neg_binomial_log_likelihood")
){
    if (is.character(error)) {
        error <- match.arg(error)
        error_fun <- paste0("get_", error)
    } else if (is.function(error)) {
        error_fun <- error
    } else {
        stop("`error` must be a character or a function.")
    }
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
        ptk_hypertune$fit_obj_list <- fit_obj_list
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
        ptk_hypertune$min_lambda <- ptk_hypertune$lambda[ptk_hypertune$lambda_min_index]
        ptk_hypertune$val_predict <- ptk_hypertune$val_predict_list[[ptk_hypertune$lambda_min_index]]
        ptk_hypertune$val_error <- error_v
        ptk_hypertune$min_error <- error_v[ptk_hypertune$lambda_min_index]
        ptk_hypertune$error_name <- error_fun
        class(ptk_hypertune) <- "ptk_hypertune"
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
    predict(object = object$fit_obj_list[[lambda_index]], newx = newx, ...)
}