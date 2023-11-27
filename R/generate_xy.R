#' @title Generate the predictor matrix in a model-specific way
#' @description Generate the numeric predictor matrix from the expression and
#' (possibly) pheno data for a certain model
#' @param expr_mat numeric matrix. The expression data, with patients as rows
#' and genes as columns.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param patient_id_col string. The name of the column in `pheno_tbl` that holds
#' the patient identifiers. Default is `"patient_id"`.
#' @param include_from_continuous_pheno vector of strings. The names of the continuous 
#' variables in `pheno_tbl to include in the predictor matrix. The values
#' will be coerced to numeric. Default is `NULL`, which means no continuous pheno
#' variables will be included.
#' @param include_from_discrete_pheno vector of strings. The names of the discrete
#' variables in `pheno_tbl` to be include in the predictor matrix. A discrete
#' variable with n levels will be converted to n-1 dummy binary variables. Default is
#' `NULL`, which means no discrete pheno variables will be included.
#' @return A numeric matrix with patients as rows and variables as columns.
#' @export
generate_predictor <- function(
    expr_mat,
    pheno_tbl,
    patient_id_col,
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL
){
    x <- expr_mat
    patient_ids <- rownames(expr_mat) # store for later
    bind_continuous <- NULL
    bind_discrete <- NULL

    check_consistent_patient_ids(
        stage = "before_generate_predictor", 
        expr = expr_mat, 
        pheno = pheno_tbl,
        patient_id_col = patient_id_col
    )

    # continuous pheno first
    if(!is.null(include_from_continuous_pheno)){
        check_tbl_columns_exist(pheno_tbl, "pheno_tbl", include_from_continuous_pheno)
        bind_continuous <- pheno_tbl[, include_from_continuous_pheno, drop = FALSE] |> 
            as.matrix()
    }
    # discrete pheno second
    if(!is.null(include_from_discrete_pheno)){
        check_tbl_columns_exist(pheno_tbl, "pheno_tbl", include_from_discrete_pheno)
        bind_discrete <- pheno_tbl[, include_from_discrete_pheno, drop = FALSE] |>
            tibble_to_binary()
    }

    # combine into numeric matrix, the predictor matrix    
    x <- x |> cbind(bind_continuous, bind_discrete)
    rownames(x) <- patient_ids

    return(x)
}


#' @title Generate the response matrix in a response-specific way
#' @description Generate the numeric response matrix from the pheno data for a
#' certain model
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param response_type string. For which response type to prepare. One of
#'  * `"survival_censored"`: response matrix (return value) will have two columns, 
#' the first is the survival time, the second is the censoring status (1 = censored, 
#' 0 = not censored)
#'  * `"binary"`: discretize the response via pfs <= `pfs_leq` (see details)
#' @param patient_id_col string. The name of the column in the pheno data file that holds 
#' patient identifiers. Default is `"patient_id"`.
#' @param pfs_col string. The name of the column in the pheno data file that holds the
#' progression-free survival (PFS) data. Default is `"pfs_years"`.
#' @param progression_col string. The name of the column in the `pheno_tbl` that holds
#' the progression status encoded as 1 = progession, 0 = no progression. 
#' Default is `"progression"`.
#' @param pfs_leq numeric. Categorize patients with progression-free survival (PFS) less than 
#' or equal `pfs_leq` as high-risk. Only used if `response_type == "binary"`. 
#' Default is 2.0.
#' @return Response matrix: a numeric matrix with patients as rows and variables as columns.
#' @details If `response_type == "binary"`, the response matrix will have one column filled 
#' with 
#' * `1` if progress is observed at a time <= `pfs_leq`,
#' * `0` if progress or censoring is observed at a time > `pfs_leq`,
#' * `NA` if censoring without progression is observed at a time <= `pfs_leq`.
#' @export
generate_response <- function(
    pheno_tbl,
    response_type,
    patient_id_col,
    pfs_col,
    progression_col,
    pfs_leq
){
    if(!(response_type %in% c("survival_censored", "binary"))){
        stop("Response type ", response_type, " is not supported")
    }
    if(response_type == "binary"){
        # remove patients consored before pfs_leq
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