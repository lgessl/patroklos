#' @title Prepare data for model fitting and predicting
#' @description Read data from files, generate the predictor and
#' respondent matirx and control quality.
#' @inherit read
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
#' @param patient_id_col string. The name of the column in the pheno data file
#' that holds the patient identifiers. Default is `"patient_id"`.
#' @param gene_id_col string. The name of the column in the expression data file
#' that holds the gene identifiers. Default is `"gene_id"`.
#' @param pfs_leq numeric. Only necessary if model discretizes response (PFS). The value 
#' of progression-free survival (PFS) above which samples are considered high-risk. 
#' Default is 2.0.
#' @return A list with two elements, `x` and `y`. `x` is the predictor matrix, `y` is
#' the response matrix with the samples as rows.
#' @export
prepare <- function(
    directory,
    model,
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    patient_id_col = "patient_id",
    gene_id_col = "gene_id",
    pfs_leq = 2.0
){
    tbls <- read(
        directory = directory,
        expr_fname = expr_fname,
        pheno_fname = pheno_fname,
        patient_id_col = patient_id_col,
        gene_id_col = gene_id_col
    )
    x <- generate_predictor(
        expr_df = tbls[["expr"]],
        pheno_df = tbls[["pheno"]],
        include_from_continuous_pheno = include_from_continuous_pheno,
        include_from_discrete_pheno = include_from_discrete_pheno
    )
    y <- generate_response(
        tbls[["pheno"]],
        model
    )
    prepare_qc(
        x = x,
        y = y
    )
    return(list("x" = x, "y" = y))
}
