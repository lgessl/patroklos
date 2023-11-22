#' @title Generate the predictor matrix in a model-specific way
#' @description Generate the numeric predictor matrix from the expression and
#' (possibly) pheno data for a certain model
#' @param expr_mat numeric matrix. The expression data, with patients as rows
#' and genes as columns.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
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


response_mapper <- list(
    "cox_lasso_zerosum" = c("pfs_yrs", "progression"),
    "lasso_zerosum" = "pfs_leq"
)

#' @title Generate the response matrix in a model-specific way
#' @description Generate the numeric response matrix from the pheno data for a
#' certain model
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
#' columns.
#' @param model string. The name of the model. One of `c("lasso-zerosum", "cox-lasso-zerosum")`.
#' @param pfs_leq numeric. Categorize patients with progression-free survival (PFS) less than 
#' or equal `pfs_leq` as high-risk. Only necessary if model discretizes response (PFS). 
#' Default is 2.0.
#' @param patient_id_col string. The name of the column in the pheno data file that holds 
#' patient identifiers. Default is `"patient_id"`.
#' @param pfs_col string. The name of the column in the pheno data file that holds the
#' progression-free survival (PFS) data. Default is `"pfs_years"`.
#' @return A numeric matrix with patients as rows and variables as columns.
#' @export
generate_response <- function(
    pheno_tbl,
    model,
    pfs_leq = 2.0,
    patient_id_col = "patient_id",
    pfs_col = "pfs_years"
){
    use <- response_mapper[[model]]
    y <- NULL
    if(length(use) == 1 && use == "pfs_leq"){ # lasso-zerosum
        # remove patients consored before pfs_leq
        na_bool <- (pheno_tbl[[pfs_col]] <= pfs_leq) & (pheno_tbl[["progression"]] == 0)
        y <- pheno_tbl[[pfs_col]] <= pfs_leq
        y <- as.numeric(y)
        dim(y) <- c(length(y), 1)
        rownames(y) <- pheno_tbl[[patient_id_col]]
        colnames(y) <- stringr::str_c("pfs_leq_", round(pfs_leq, 1))
        y[na_bool, ] <- NA
    } else { # cox-lasso-zerosum
        y <- pheno_tbl[, use] |> as.matrix()
        rownames(y) <- pheno_tbl[[patient_id_col]]
    }
    return(y)
}