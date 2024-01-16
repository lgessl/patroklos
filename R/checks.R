#' @title Check if a tibble has certain columns
#' @param tbl tibble. The tibble to check.
#' @param tbl_name string. How to refer to the tibble in an error message.
#' @param col_names vector of strings. The names of the columns to check for.
#' @return NULL. Throws an error if one of the columns is not in the tibble.
check_tbl_columns_exist <- function(
    tbl,
    tbl_name,
    col_names
){
    for(cname in col_names){
        if(is.null(tbl[[cname]])){
            stop("Column ", cname, " does not exist in ", tbl_name, ".")
        }
    }
}

check_consistent_patient_ids <- function(
    stage = c("preprocessing", "before_generate_predictor", "after_generate_xy"),
    expr,
    pheno,
    data_spec = NULL
){
    # Extract
    stage <- match.arg(stage)
    patient_id_col <- data_spec$patient_id_col
    patient_ids_expr <- NULL
    patient_ids_pheno <- NULL
    # case 1: during preprocessing
    if(stage == "preprocessing"){
        if(!is.data.frame(expr)) stop("`expr` must be a tibble")
        if(!is.data.frame(pheno)) stop("`pheno` must be a tibble")
        patient_ids_expr <- colnames(expr)[-1]
        patient_ids_pheno <- pheno[[patient_id_col]]
    }
    # case 2: right before generate_predictor
    if(stage == "before_generate_predictor"){
        if(!is.matrix(expr)) stop("`expr` must be a matrix")
        if(!is.data.frame(pheno)) stop("`pheno` must be a tibble")
        patient_ids_expr <- rownames(expr)
        patient_ids_pheno <- pheno[[patient_id_col]]
    }
    # case 3: after generate_predictor and generate_response
    if(stage == "after_generate_xy"){
        if(!is.matrix(expr)) stop("`expr` must be a matrix")
        if(!is.matrix(pheno)) stop("`pheno` must be a matrix")
        patient_ids_expr <- rownames(expr)
        patient_ids_pheno <- rownames(pheno)
    }

    if(is.null(patient_ids_expr))
        stop("Expression data does not have sample ids")
    if(is.null(patient_ids_pheno))
        stop("Pheno data does not have sample ids")
    if(!identical(patient_ids_expr, patient_ids_pheno)){
        stop("Patient ids in expression and pheno data are not identical.")
    }
}


check_available <- function(
    x,
    y
){
    if(any(is.na(x)) || any(is.na(y))){
        stop("There are missing values in the predictor or response matrix.")
    }
}


check_fitter <- function(
    fitter,
    optional_args = NULL
){
    if(!is.function(fitter)){
        stop("fitter must be a function.")
    }
    if(!(is.null(optional_args) || is.list(optional_args))){
        stop("optional_args must be a list.")
    }
    parameters <- names(formals(fitter))
    for(pos_arg in c("x", "y")){
        if(!(pos_arg %in% parameters)){
            stop("fitter must have an argument named '", pos_arg, "'.")
        }
    }
    if(!("..." %in% parameters)){
        for(arg in optional_args){
            if(!(arg %in% parameters)){
                stop("fitter must have an argument named '", arg, "'.")
            }
        }
    }
}