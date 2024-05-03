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

ptk_ranger <- function(x, y, ...){
    args <- list(...)
    args <- args[names(args) != "li_var_suffix"]
    ptk_ranger_obj <- do.call(
        ranger::ranger, 
        c(list("x" = x, "y" = y), args)
    )
    class(ptk_ranger_obj) <- c("ptk_ranger", class(ptk_ranger_obj))
}

predict.ptk_ranger <- function(object, newx, ...){
    y <- do.call(
        ranger::predict.ranger, 
        c(list("object" = object, "data" = newx), list(...))
    )[[1]]
    dim(y) <- NULL
    names(y) <- rownames(newx)
    y
}