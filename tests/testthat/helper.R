# Generate random mock data
generate_mock_data <- function(
    directory = "data",
    n_samples = 10,
    n_genes = 5,
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    patient_id_col = "patient_id",
    gene_id_col = "gene_id"
){
    if(file.exists(file.path(directory, expr_fname)) &&
        file.exists(file.path(directory, pheno_fname))){
    }
    # expression
    expr_mat <- matrix(
        sample(1:100, n_samples*n_genes, replace = TRUE),
        nrow = n_genes
    )
    colnames(expr_mat) <- stringr::str_c("sample_", 1:n_samples)
    expr_tbl <- tibble::as_tibble(expr_mat)
    expr_tbl[[gene_id_col]] <- stringr::str_c("gene_", 1:n_genes)
    # pheno
    pheno <- tibble::tibble(.rows = n_samples)
    pheno[["progression"]] <- sample(
        c(0, 1),
        size = n_samples,
        replace = TRUE
    )
    pheno[["pfs_yrs"]] <- round(rnorm(n_samples, 2, 1), 2)
    pheno[[patient_id_col]] <- colnames(expr_mat)
    # write to csv
    readr::write_csv(
        expr_tbl,
        file.path(directory, expr_fname)
    )
    readr::write_csv(
        pheno,
        file.path(directory, pheno_fname)
    )
}