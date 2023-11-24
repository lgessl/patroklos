#' @title Prepare data for model fitting and predicting
#' @description Provided expression matrix and pheno tibble, generate the 
#' predictor and respondent matirx model-specifically and check.
#' @param expr_mat numeric matrix. The expression data, with patients as rows
#' and genes as columns. The row names must be unique patient identifiers.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
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
#' @param patient_id_col string. The name of the column in `pheno_tbl`
#' that holds the patient identifiers. Default is `"patient_id"`.
#' @param pfs_col string. The name of the column in `pheno_tbl` that holds the
#' progression-free survival (PFS) values. Default is `"pfs_years"`.
#' @param progression_col string. The name of the column in `pheno_tbl` that
#' holds the progression status encoded as 1 = progression, 0 = no progression. 
#' Default is `"progression"`.
#' @param pfs_leq numeric. Only necessary if model discretizes response (PFS). The value 
#' of progression-free survival (PFS) above which samples are considered high-risk. 
#' Default is 2.0.
#' @return A list with two numeric matrices, `x` and `y`. `x` is the predictor matrix, 
#' `y` is the response matrix, both ready to for model fitting or predicting.
#' @export
prepare <- function(
    expr_mat,
    pheno_tbl,
    model,
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    patient_id_col = "patient_id",
    pfs_col = "pfs_years",
    progression_col = "progression",
    pfs_leq = 2.0
){
    x <- generate_predictor(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        patient_id_col = patient_id_col,
        include_from_continuous_pheno = include_from_continuous_pheno,
        include_from_discrete_pheno = include_from_discrete_pheno
    )

    y <- generate_response(
        pheno_tbl = pheno_tbl,
        model = model,
        patient_id_col = patient_id_col,
        pfs_col = pfs_col,
        progression_col = progression_col,
        pfs_leq = pfs_leq
    )

    x_y <- ensure_available(
        x = x,
        y = y
    )
    x <- x_y[["x"]]
    y <- x_y[["y"]]

    qc_prepare(
        x = x,
        y = y
    )
    return(list("x" = x, "y" = y))
}
