# Some straightforward wrappers around existing models

#' @title Call [`zeroSum::zeroSum()`] with zeroSum.weights corresponding to features 
#' from pheno data set to zero
#' @description Recognize features added to predcitor from pheno via a common suffix
#' and set zeroSum.weights to 0 for those features (since they are, unlike RNAseq data, 
#' not affected by rescaling), then call [`zeroSum::zeroSum()`]
#' @param x See [`zeroSum::zeroSum()`].
#' @param y See [`zeroSum::zeroSum()`].
#' @param pheno_regexp string that matches exactly those features (column names) in 
#' `x` that were added from pheno data. Default is `"\\+\\+$"`.
#' @param ... further arguments passed to [`zeroSum::zeroSum()`].
#' @return `zeroSumFit` S3 object as returned by [`zeroSum::zeroSum()`].
#' @export
zeroSumEI <- function(
    x,
    y,
    pheno_regexp = "\\+\\+$",
    ...
){
    zeroSum.weights <- as.numeric(!stringr::str_detect(colnames(x), pheno_regexp))
    fit_obj <- zeroSum::zeroSum(
        x = x, 
        y = y, 
        zeroSum.weights = zeroSum.weights,
        ...
        )
    fit_obj$zeroSum.weights <- zeroSum.weights
    return(fit_obj)
}

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