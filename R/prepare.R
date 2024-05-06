data_prepare <- function(
    self,
    private,    
    model
){
    if(length(model$split_index) != 1){
        stop("Model must specify exactly one split index, but we have ", 
            length(model$split_index))
    }
    if(length(model$time_cutoffs) != 1){
        stop("Model must specify exactly one time cutoff")
    }
    if(is.null(self$cohort)){
        stop("Cohort must be set to `'train'` or `'test'`")
    }

    # Subset data to cohort
    split_colname <- paste0(self$split_col_prefix, model$split_index)
    patient_id_col <- self$patient_id_col
    if(!split_colname %in% colnames(self$pheno_tbl))
        stop("Column ", split_colname, " not found in pheno table.")
    cohort_bool <- self$pheno_tbl[[split_colname]] == self$cohort
    if(all(cohort_bool) || all(!cohort_bool))
        stop("All patients are in the same cohort")
    pheno_tbl <- self$pheno_tbl[cohort_bool, ]
    expr_mat <- self$expr_mat[pheno_tbl[[patient_id_col]], , drop = FALSE]

    x <- generate_predictor(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data = self,
        model = model
    )
    y <- generate_response(
        pheno_tbl = pheno_tbl,
        data = self,
        model = model
    )

    return(list("x" = x, "y" = y))
}
