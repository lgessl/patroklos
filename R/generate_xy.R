#' @title Generate the predictor matrix in a model-specific way
#' @description Generate the numeric predictor matrix from the expression and
#' (possibly) pheno data for a certain model
#' @param expr_mat numeric matrix. The expression data, with patients as rows
#' and genes as columns.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param data_spec DataSpec S3 object. Specifications on the data (`expr_mat` and 
#' `model_spec`). See the the constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model to prepare for. 
#' See the the constructor `ModelSpec()` for details.
#' @return A numeric matrix with patients as rows and variables as columns.
#' @export
generate_predictor <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec
){
    # Extract
    patient_id_col <- data_spec$patient_id_col
    include_from_continuous_pheno <- model_spec$include_from_continuous_pheno
    include_from_discrete_pheno <- model_spec$include_from_discrete_pheno

    x <- expr_mat
    patient_ids <- rownames(expr_mat) # store for later
    bind_continuous <- NULL
    bind_discrete <- NULL

    check_consistent_patient_ids(
        stage = "before_generate_predictor", 
        expr = expr_mat, 
        pheno = pheno_tbl,
        data_spec = data_spec
    )

    # continuous pheno first
    if(!is.null(include_from_continuous_pheno)){
        check_tbl_columns_exist(pheno_tbl, "pheno_tbl", include_from_continuous_pheno)
        bind_continuous <- pheno_tbl[, include_from_continuous_pheno, drop = FALSE] |> 
            as.matrix()
        colnames(bind_continuous) <- colnames(bind_continuous) |> 
            stringr::str_c(model_spec$append_to_includes)
    }
    # discrete pheno second
    if(!is.null(include_from_discrete_pheno)){
        check_tbl_columns_exist(pheno_tbl, "pheno_tbl", include_from_discrete_pheno)
        bind_discrete <- pheno_tbl[, include_from_discrete_pheno, drop = FALSE] |>
            tibble_to_binary()
        colnames(bind_discrete) <- colnames(bind_discrete) |>
            stringr::str_c(model_spec$append_to_includes)
    }

    # combine into numeric matrix, the predictor matrix    
    x <- x |> cbind(bind_continuous, bind_discrete)
    rownames(x) <- patient_ids

    return(x)
}


#' @title Generate the response matrix in a model-specific way
#' @description Generate the numeric response matrix from the pheno data for a
#' certain model
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param data_spec DataSpec S3 object. Specifications on the data (`pheno_tbl`). See
#' the constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model to prepare for.
#' See the constructor `ModelSpec()` for details.
#' @return Response matrix: a numeric matrix with patients as rows and variables as columns.
#' @details If `model_spec$response_type == "binary"`, the response matrix will have one 
#' column filled with 
#' * `1` if progress is observed at a time <= `pfs_leq`,
#' * `0` if progress or censoring is observed at a time > `pfs_leq`,
#' * `NA` if censoring without progression is observed at a time <= `pfs_leq`.
#' @export
generate_response <- function(
    pheno_tbl,
    data_spec,
    model_spec
){

    # Extract
    pfs_col <- data_spec$pfs_col
    progression_col <- data_spec$progression_col
    patient_id_col <- data_spec$patient_id_col
    pfs_leq <- model_spec$pfs_leq
    response_type <- model_spec$response_type

    if(!(response_type %in% c("survival_censored", "binary"))){
        stop("Response type ", response_type, " is not supported")
    }

    if(response_type == "binary"){
        # flag patients consored before pfs_leq as NA
        na_bool <- (pheno_tbl[[pfs_col]] <= pfs_leq) & (pheno_tbl[[progression_col]] == 0)
        y <- pheno_tbl[[pfs_col]] <= pfs_leq
        y <- as.numeric(y)
        dim(y) <- c(length(y), 1)
        rownames(y) <- pheno_tbl[[patient_id_col]]
        colnames(y) <- stringr::str_c("pfs_leq_", round(pfs_leq, 1))
        y[na_bool, ] <- NA
    } else if(response_type == "survival_censored"){
        y <- pheno_tbl[, c(pfs_col, progression_col)] |> as.matrix()
        rownames(y) <- pheno_tbl[[patient_id_col]]
    }

    return(y)
}