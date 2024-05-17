#' @title Wrap [`ranger::ranger()`] into a patroklos-compliant fitter 
#' @description This function is a patroklos-compliant fitter with validated 
#' predictions, it has a return value with a `oob_predict` attribute.
#' @inheritParams ranger::ranger
#' @param skip_on_invalid_input Logical. If `TRUE` and invalid input is detected,
#' return `"next"` instead of an error. This is useful when calling this function 
#' from inside [`nested_pseudo_cv()`].
#' @param ... Further arguments passed to the wrapped function.
#' @return A `ptk_ranger` S3 object, a `ranger` S3 object with the `predictions` 
#' attribute renamed to `oob_predict`.
#' @export
ptk_ranger <- function(x, y, mtry, skip_on_invalid_input = FALSE, ...){
    if (ncol(x) < mtry) {
        if (skip_on_invalid_input) return("next")
        stop("mtry must be less than the number of features.")
    }
    ptk_ranger_obj <- ranger::ranger(x = x, y = y, mtry = mtry, ...)
    # Rename OOB predictions
    ptk_ranger_obj$oob_predict <- ptk_ranger_obj$predictions
    ptk_ranger_obj$predictions <- NULL
    class(ptk_ranger_obj) <- c("ptk_ranger", class(ptk_ranger_obj))
    ptk_ranger_obj
}

#' @title Wrap [`ranger::predict.ranger()`] into a patroklos-compliant predict
#' function
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
#' @description This function is a patroklos-compliant fitter with integrated 
#' CV and, if `length(lambda) == 1`, also a patroklos-compliant fitter with 
#' validated predictions.
#' @inheritParams zeroSum::zeroSum
#' @param exclude_pheno_from_lasso Logical. If `TRUE`, set LASSO penalty weights
#' corresponding to features from the pheno data to zero.
#' @param binarize_predictions numeric or NULL. If not NULL, the predict method 
#' for the returned `ptk_zerosum` object will binarize the predictions using the
#' `binarize_predictions` as a threshold.
#' @return A `ptk_zerosum` S3 object. 
#' @details A *patroklos-compliant fitter with integrated CV* is a patroklos-
#' compliant fitter (cf. details of [`ptk_ranger()`]) whose return value has an 
#' attribute `cv_predict_list`, a list of numeric vectors holding cross-validated 
#' predictions for every value of the model hyperparameter lambda.
#' @export
ptk_zerosum <- function(
    x,
    y,
    exclude_pheno_from_lasso = TRUE,
    binarize_predictions = NULL,    
    ...,
    nFold = 10,
    zeroSum.weights = NULL, # Do not partially match these
    penalty.factor = NULL,
    family = "binomial"
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
    fit_obj <- zeroSum::zeroSum(x = x, y = y, nFold = nFold, 
        zeroSum.weights = zeroSum.weights, penalty.factor = penalty.factor, 
        family = family, ...)
    fit_obj$cv_predict_list <- fit_obj$cv_predict
    fit_obj$cv_predict <- NULL
    if (family == "binomial") {
        fit_obj$cv_predict_list <- lapply(fit_obj$cv_predict_list,
            function(v) 1/(1+exp(-v)))
    }
    if (is.numeric(binarize_predictions))
        fit_obj$cv_predict_list <- lapply(fit_obj$cv_predict_list,
            function(v) as.numeric(v > binarize_predictions))
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
#' @return A numeric vector with the prediction for every sample. Unline 
#' `zeroSum::predict.zeroSum()`, 
#' * if `object$type == 2` (i.e., `object` was fitted with `family = "binomial"
#' by `ptk_zerosum`), the predictions are binomial probabilities,
#' * if `object$type == 4` (i.e., `object` was fitted with `family = "cox" by
#' `ptk_zerosum`), the predictions are the exponential of the linear predictor, 
#' i.e., `exp(beta_0 + beta^T x)`.
#' * If `object$binarizePredictions` is not `NULL`, the predictions are binarized
#' via `as.numeric(y > object$binarizePredictions)`.
#' @export
predict.ptk_zerosum <- function(
    object,
    newx,
    ...
){
    class(object) <- "zeroSum"
    args <- list(...)
    if(!is.null(args[["type"]]))
        stop("`type` argument is not allowed.") 
    type <- NULL
    if(object$type %in% c(2, 4)) # 2 \mapsto binomial, 4 \mapsto cox
        # get probabilities (binomial), exp(beta_0 + beta^T x) (cox)
        type <- "response"     
    y <- predict(object = object, newx = newx, type = type, ...)
    if (!is.null(object$binarizePredictions))
        y <- as.numeric(y > object$binarizePredictions)
    return(y)
}