#' @title Prepare data for model fitting and predicting
#' @description Read data from files, generate the predictor and
#' respondent matirx and control quality.
#' @param directory string. The directory where both expression and 
#'  pheno data files lie.
#' @param model string. For which model to prepare. One of
#'  * `"cox_lasso_zerosum"` (Cox proportional hazards with LASSO reularization 
#'  and zero-sum constraint)
#'  * `"lasso_zerosum"` (ordinary linear regression with LASSO regularization)
#' @param include_from_continuous_pheno vector of strings. The names of the 
#'  *continuous* variables in the pheno data file to be included in the
#'  predictor matrix. The values will be coerced to numeric. Default is `NULL`,
#'  which means no continuous pheno variables will be included.
#' @param include_from_discrete_pheno vector of strings. The names of the
#'  *discrete* variables in the pheno data file to be included in the predictor
#'  matrix. A discrete variable with n levels will be converted to n-1 binary
#'  variables. Default is `NULL`, which means no discrete pheno variables will
#'  be included.
#' @param expr_fname string. The name of the expression .csv file inside `directory`.
#'  Default is `"expr.csv"`.
#' @param pheno_fname string. The name of the pheno data .csv inside `directory`.
#'  Default is `"pheno.csv"`.
#' @details The pheno .csv files holds the samples as rows (with the unique sample names
#'  in the very first column), the variables as columns. We need at least the columns
#'  * progression (integer, 0 for censored, 1 for uncensored),
#'  * pfs_yrs (integer, the time in years to progression or censoring),
#'  * `ìnclude_from_continuous_pheno`,
#'  * `include_from_discrete_pheno`.
#' The expr .csv file holds the genes as rows (with the unique gene names in the very
#' first column), the samples as columns. The sample names in pheno and expression files
#' must be identical and in the same order.
#' @return A list with two elements, `x` and `y`. `x` is the predictor matrix, `y` is
#' the respondent matrix with the samples as rows.
prepare <- function(
    directory,
    model,
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv"
){
    tbls <- read(
        directory = directory,
        model = model,
        include_from_continuous_pheno = include_from_continuous_pheno,
        include_from_discrete_pheno = include_from_discrete_pheno
    )
    x <- generate_predictor(
        expr_df = tbls[["expr"]],
        pheno_df = tbls[["pheno"]],
        include_from_continuous_pheno = include_from_continuous_pheno,
        include_from_discrete_pheno = include_from_discrete_pheno
    )
    prepare_model <- switch(model,
        "cox_lasso_zerosum" = prepare_cox_lasso_zerosum,
        "lasso_zerosum" = prepare_cox_lasso)
    y <- generate_respondent(
        tbls[["pheno"]]
    )
    prepare_qc(
        x = x,
        y = y
    )
    return(list("x" = x, "y" = y))
}
