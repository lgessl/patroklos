#' @title Build a fitter with integrated cross-validation from a fitter
#' @description This function operator takes a patroklos-compliant fitter and 
#' builds a patroklos-compliant fitter with integrated cross validation from 
#' it.
#' @param fitter A patroklos-compliant fitter.
#' @param metric A character string or a function. If a character string, it
#' must be one of "accuracy", "roc_auc", or "binomial_log_likelihood". If a
#' function, it must take two arguments: `y` and `y_hat`, and return a numeric
#' scalar. The returned fitter will use it to calculate the goodness of validated 
#' predictions.
#' @return A patroklos-compliant fitter with integrated cross-validation.
integrated_cv <- function(
    fitter, 
    metric = c("accuracy", "roc_auc", "binomial_log_likelihood")
){
    if (is.character(metric)) {
        metric <- match.arg(metric)
        metric_fun <- paste0("get_", metric)
    } else if (is.function(metric)) {
        metric_fun <- metric
    } else {
        stop("`metric` must be a character or a function.")
    }
    force(fitter)
    force(metric_fun)

    cv_fitter <- function(x, y, ...) {
        # Get all possible combinations of hyperparameters
        grid <- expand.grid(list(...), stringsAsFactors = FALSE)
        fit_obj_list <- lapply(seq(nrow(grid)), function(i) {
            do.call(fitter, c(list(x = x, y = y), as.list(grid[i, ])))
        })
        # Build ptk_cv S3 object
        ptk_cv <- list()
        ptk_cv$fit_obj_list <- fit_obj_list
        # Promote validated predictions to top level
        ptk_cv$val_predict_list <- lapply(fit_obj_list, function(x) x$val_predict)
        ptk_cv$lambda <- apply(grid, 1, function(r) paste(names(grid), r, 
            sep = "=", collapse = ", "))
        # Evaluate according to metric
        metric_v <- sapply(ptk_cv$val_predict_list, function(y_hat) 
            do.call(metric_fun, list(y, y_hat)))
        ptk_cv$best_lambda_index <- which.max(metric_v)
        ptk_cv$best_lambda <- ptk_cv$lambda[ptk_cv$best_lambda_index]
        ptk_cv$val_predict <- ptk_cv$val_predict_list[[ptk_cv$best_lambda_index]]
        ptk_cv$val_metric <- metric_v
        ptk_cv$best_metric <- metric_v[ptk_cv$best_lambda_index]
        ptk_cv$metric_name <- metric_fun
        class(ptk_cv) <- "ptk_cv"
        ptk_cv
    }
}

#' @export
predict.ptk_cv <- function(
    object,
    newx,
    lambda_index = object$best_lambda_index,
    ...
){
    predict(object = object$fit_obj_list[[lambda_index]], newx = newx, ...)
}