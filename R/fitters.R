# Some straightforward wrappers around existing models

#' @title Call [zeroSum::zeroSum()] with zeroSum.weights corresponding to features 
#' from pheno data set to zero
#' @description Recognize features added to predcitor from pheno via a common suffix
#' and set zeroSum.weights to 0 for those features (since they are, unlike RNAseq data, 
#' not affected by rescaling), then call [zeroSum::zeroSum()]
#' @param x See [zeroSum::zeroSum()].
#' @param y See [zeroSum::zeroSum()].
#' @param pheno_regexp string that matches exactly those features (column names) in 
#' `x` that were added from pheno data. Default is `"\\+\\+$"`.
#' @param ... further arguments passed to [zeroSum::zeroSum()].
#' @return `zeroSumFit` S3 object as returned by [zeroSum::zeroSum()].
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