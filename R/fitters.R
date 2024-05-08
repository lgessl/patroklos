#' @title Wrap [`ranger::ranger()`] into a patroklos-compliant fit function
#' @description You can now use it as the `fitter` attribute of a `Model` object.
#' @param x Predictor data as a named numeric matrix. Samples correspond to rows,
#' an attribute `li_var_suffix` is expected to be present.
#' @param y Response data as numeric vector. The length of `y` must match the
#' number of rows in `x`.
#' @param ... Further arguments passed to the wrapped function.
#' @return A `ptk_ranger` S3 object.
#' @details The interface of this function is the interface of a patroklos-compliant
#' fitter.
#' @export
ptk_ranger <- function(x, y, ...){
    ptk_ranger_obj <- ranger::ranger(x = x, y = y, ...)
    class(ptk_ranger_obj) <- c("ptk_ranger", class(ptk_ranger_obj))
    ptk_ranger_obj
}

#' @title Wrap [`ranger::predict.ranger()`] into a patroklos-compliant predict
#' function
#' @description `Model$predict()` uses this function as a predict method for
#' the `predict()` generic.
#' @param object A fit S3 object.
#' @param newx Predictor data as a named numeric matrix. Samples correspond to 
#' rows, an attribute `li_var_suffix` is expected to be present.
#' @param ... Further arguments passed to the wrapped function.
#' @return A named numeric vector of predictions. 
#' @details The interface of this function is the interface of a patroklos-compliant
#' predict method.
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
#' @inheritParams ptk_ranger
#' @param exclude_pheno_from_lasso Logical. If `TRUE`, set LASSO penalty weights
#' corresponding to features from the pheno data to zero.
#' @param binarize_predictions numeric or NULL. If not NULL, the predict method 
#' for the returned `ptk_zeroSum` object will binarize the predictions using the
#' `binarize_predictions` as a threshold.
#' @return A `ptk_zeroSum` S3 object. 
#' @export
ptk_zeroSum <- function(
    x,
    y,
    exclude_pheno_from_lasso = TRUE,
    binarize_predictions = NULL,    
    ...
){
    early_bool <- get_early_bool(x)
    zerosum_weights <- as.numeric(early_bool) 
    penalty_factor <- rep(1, ncol(x))
    if (exclude_pheno_from_lasso)
        penalty_factor <- as.numeric(early_bool)
    fit_obj <- zeroSum::zeroSum(
        x = x, 
        y = y, 
        zeroSum.weights = zerosum_weights,
        penalty.factor = penalty_factor,
        ...
        )
    fit_obj$zeroSum.weights <- zeroSum.weights
    fit_obj$binarize_predictions <- binarize_predictions
    class(fit_obj)[1] <- "ptk_zeroSum"
    return(fit_obj)
}