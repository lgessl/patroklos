#' @title Prepare data for model fitting and predicting
#' @description Provided expression matrix and pheno tibble, generate the 
#' predictor and response matrix response-type-dependently and do checks. At
#' this point, the ModelSpec must specifiy exactly one time cutoff and exactly
#' one split index.
#' @param expr_mat numeric matrix. The expression data, with patients as rows
#' and genes as columns. The row names must be unique patient identifiers.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param data_spec DataSpec S3 object. Specifications on the data (`expr_mat`, 
#' `pheno_tbl`). See the the constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model to prepare for. 
#' See the the constructor `ModelSpec()` for details.
#' @return A list with two numeric matrices, `x` and `y`. `x` is the predictor matrix, 
#' `y` is the response matrix, both ready to for model fitting or predicting. Row names 
#' match.
#' @export
prepare <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec
){
    if(length(model_spec$split_index) != 1){
        stop("ModelSpec must specify exactly one split index")
    }
    if(length(model_spec$time_cutoffs) != 1){
        stop("ModelSpec must specify exactly one time cutoff")
    }
    if(is.null(data_spec$cohort)){
        stop("data_spec$cohort must be set to `'train'` or `'test'`")
    }

    # Subset data to cohort
    split_colname <- paste0(data_spec$split_col_prefix, model_spec$split_index)
    patient_id_col <- data_spec$patient_id_col
    if(!split_colname %in% colnames(pheno_tbl))
        stop("Column ", split_colname, " not found in pheno table.")
    cohort_bool <- pheno_tbl[[split_colname]] == data_spec$cohort
    if(all(cohort_bool) || all(!cohort_bool))
        stop("All patients are in the same cohort")
    pheno_tbl <- pheno_tbl[cohort_bool, ]
    expr_mat <- expr_mat[pheno_tbl[[patient_id_col]], , drop = FALSE]

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
    
    return(list("x" = x, "y" = y))
}
