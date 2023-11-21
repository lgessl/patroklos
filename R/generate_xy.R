mapper <- list(
    "cox_lasso_zerosum" = c("pfs_yrs", "progression"),
    "lasso_zerosum" = "pfs_leq"
)

generate_predictor <- function(
    expr_mat,
    pheno_tbl,
    include_from_continuous_pheno,
    include_from_discrete_pheno,
    gene_id_col = "gene_id"
){
    x <- expr_mat
    patient_ids <- rownames(expr_mat)
    bind_continuous <- pheno_tbl[, include_from_continuous_pheno, drop = FALSE] |> 
        as.matrix()
    bind_discrete <- pheno_tbl[, include_from_discrete_pheno, drop = FALSE] |>
        tibble_to_binary()
    x <- x |> cbind(bind_continuous, bind_discrete)
    rownames(x) <- patient_ids
    return(x)
}

# generate the response matrix or vector in a model-specific way
generate_response <- function(
    pheno_tbl,
    model,
    pfs_leq = 2.0,
    patient_id_col = "patient_id"
){
    use <- mapper[[model]]
    y <- NULL
    if(length(use) == 1 && use == "pfs_leq"){
        # remove patients consored before pfs_leq
        rm_bool <- (pheno_tbl[["pfs_yrs"]] <= pfs_leq) & (pheno_tbl[["progression"]] == 0)
        y <- pheno_tbl[["pfs_yrs"]] <= pfs_leq
        names(y) <- pheno_tbl[[patient_id_col]]
        y <- y[!rm_bool]
    } else {
        y <- pheno_tbl[, use]
    }
    return(y)
}
