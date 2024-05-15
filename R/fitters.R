#' @title Wrap [`ranger::ranger()`] into a patroklos-compliant fit function
#' @description You can now use it as the `fitter` attribute of a `Model` object.
#' @param x Predictor data as a named numeric matrix. Samples correspond to rows,
#' an attribute `li_var_suffix` is expected to be present.
#' @param y Response data as numeric vector. The length of `y` must match the
#' number of rows in `x`.
#' @param ... Further arguments passed to the wrapped function.
#' @return A `ptk_ranger` S3 object, a `ranger` S3 object with the `predictions` 
#' attribute renamed to `oob_predict`.
#' @details A *patroklos-compliant fitter* is a function with the same parameters
#' as this function and a return value that is an S3 object with a *patroklos-
#' compliant predict method* (see details of [`predict.ptk_ranger()`]).
#' 
#' A *patroklos compliant fitter with validated predictions* is a patroklos-
#' compliant fitter whose return value has an attribute `cv_predict` or 
#' `oob_predict`, a numeric vector holding cross-validated or out-of-bag (OOB) 
#' predictions, respectively.
#' @export
ptk_ranger <- function(x, y, ...){
    ptk_ranger_obj <- ranger::ranger(x = x, y = y, ...)
    # Rename OOB predictions
    ptk_ranger_obj$oob_predict <- ptk_ranger_obj$predictions
    ptk_ranger_obj$predictions <- NULL
    class(ptk_ranger_obj) <- c("ptk_ranger", class(ptk_ranger_obj))
    ptk_ranger_obj
}

#' @title Wrap [`ranger::predict.ranger()`] into a patroklos-compliant predict
#' function
#' @description `Model$predict()` uses this function as a predict method for
#' the `predict()` generic.
#' @param object A ptk_ranger S3 object.
#' @param newx Predictor data as a named numeric matrix. Samples correspond to 
#' rows, an attribute `li_var_suffix` is expected to be present.
#' @param ... Further arguments passed to the wrapped function.
#' @return A named numeric vector of predictions from `newx`. 
#' @details A *patroklos-compliant predict method* is a function with the same
#' signature as this function. 
#' @export
predict.ptk_ranger <- function(object, newx, ...){
    class(object) <- "ranger"
    y <- predict(object, data = newx, ...)[[1]]
    dim(y) <- NULL
    names(y) <- rownames(newx)
    y
}

#' @title Wrap [`zeroSum::zeroSum()`] into a patroklos-compliant fit function
#' @description You can now use it as the `fitter` attribute of a `Model` object.
#' This fitter is suited for early integration.
#' @inheritParams zeroSum::zeroSum
#' @param exclude_pheno_from_lasso Logical. If `TRUE`, set LASSO penalty weights
#' corresponding to features from the pheno data to zero.
#' @param binarize_predictions numeric or NULL. If not NULL, the predict method 
#' for the returned `ptk_zerosum` object will binarize the predictions using the
#' `binarize_predictions` as a threshold.
#' @return A `ptk_zerosum` S3 object. 
#' @details A *patroklos-compliant fitter with integrated CV* is a patroklos-
#' compliant fitter (cf. details of [`ptk_ranger()`] whose return value has an 
#' attribute `cv_predict_list`, a list of numeric vectors holding cross-validated 
#' predictions for every value of the model hyperparameter lambda.
#' @export
ptk_zerosum <- function(
    x,
    y,
    exclude_pheno_from_lasso = TRUE,
    binarize_predictions = NULL,    
    ...,
    zeroSum.weights = NULL, # Do not partially match these
    penalty.factor = NULL
){
    stopifnot(is.logical(exclude_pheno_from_lasso))
    stopifnot(is.null(binarize_predictions) || 
        (is.numeric(binarize_predictions) && binarize_predictions > 0 && 
        binarize_predictions < 1))
    early_bool <- get_early_bool(x, for_li = FALSE)
    # zeroSum weights
    if (is.null(zeroSum.weights))
        zeroSum.weights <- as.numeric(early_bool) 
    # LASSO penalty factors
    if (exclude_pheno_from_lasso) {
        if(is.null(penalty.factor)) {
            penalty.factor <- as.numeric(early_bool)
        } else {
            penalty.factor[!early_bool] <- 0
        }
    }
    fit_obj <- zeroSum::zeroSum(x = x, y = y, zeroSum.weights = zeroSum.weights, 
        penalty.factor = penalty.factor, ...)
    fit_obj$cv_predict_list <- fit_obj$cv_predict
    fit_obj$cv_predict <- NULL
    if (length(fit_obj$cv_predict_list) == 1)
        fit_obj$cv_predict <- fit_obj$cv_predict_list[[1]]
    if (fit_obj$useZeroSum)
        fit_obj$zeroSumWeights <- zeroSum.weights
    if (fit_obj$standardize == 0)
        fit_obj$penaltyFactor <- penalty.factor
    fit_obj$binarizePredictions <- binarize_predictions
    class(fit_obj) <- "ptk_zerosum"
    return(fit_obj)
}

#' @title Wrap [`zeroSum::predict.zeroSum()`] into a patroklos-compliant predict function
#' @description `Model$predict()` uses this function as a predict method for
#' the `predict()` generic.
#' @inheritParams predict.ptk_ranger
#' @param object A ptk_patroklos S3 object.
#' @return A numeric vector with the prediction for every sample. If the 
#' `binarize_predictions` attribute of the `ptk_zerosum` object is NULL, return 
#' the linear predictor beta_0 + beta^T x. Otherwise, return the binary class 
#' labels by thresholding the logistic transformation of the linear predictor 
#' via `binarize_predictions`.
#' @export
predict.ptk_zerosum <- function(
    object,
    newx,
    ...
){
    class(object) <- "zeroSum"
    args <- list(...)
    if(!is.null(args[["type"]]))
        stop("`type` argument is not allowed. The `binarizePredictions` attribute 
            of the `ptk_zerosum` object determines the prediction type.")
    # type "link" returns linear predictor beta_0 + beta^T x
    y <- do.call("predict", c(
        list(object = object, newx = newx, type = "link"), 
        args
    ))
    if (!is.null(object$binarizePredictions))
        y <- as.numeric(1/(1+exp(-y)) > object$binarizePredictions)
    return(y)
}