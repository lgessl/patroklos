#' @title Prepare data for model fitting and predicting
#' @description Provided expression matrix and pheno tibble, generate the 
#' predictor and response matrix response-type-dependently and do checks.
#' @param expr_mat numeric matrix. The expression data, with patients as rows
#' and genes as columns. The row names must be unique patient identifiers.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param data_spec DataSpec S3 object. Specifications on the data (`expr_mat`, 
#' `pheno_tbl`). See the the constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model to prepare for. 
#' See the the constructor `ModelSpec()` for details.
#' @return A list with two numeric matrices, `x` and `y`. `x` is the predictor matrix, 
#' `y` is the response matrix, both ready to for model fitting or predicting.
#' @export
prepare <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec
){
    x <- generate_predictor(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec
    )

    y <- generate_response(
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec
    )
    check_consistent_patient_ids(
        stage = "after_generate_xy",
        expr = x,
        pheno = y
    )
    x_y <- intersect_by_names(x, y, rm_na = TRUE)
    x <- x_y[[1]]
    y <- x_y[[2]]

    qc_prepare(
        x = x,
        y = y
    )
    return(list("x" = x, "y" = y))
}
