# Generate random mock data
generate_mock_files <- function(
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


generate_mock_data <- function(
    n_samples = 10,
    n_genes = 5,
    n_na_in_pheno = 3,
    to_csv = NULL
){
    # expression matrix
    expr_mat <- matrix(
        sample(1:100, n_samples*n_genes, replace = TRUE),
        nrow = n_samples
    ) |> log()
    rownames(expr_mat) <- stringr::str_c("sample_", 1:n_samples)
    colnames(expr_mat) <- stringr::str_c("gene_", 1:n_genes)
    # to expression tibble
    expr_tbl <- expr_mat |> t() |> tibble::as_tibble(rownames = "gene_id")

    # pheno (n_samples x 5 tibble)
    pheno_tbl <- tibble::tibble(.rows = n_samples)
    pheno_tbl[["patient_id"]] <- rownames(expr_mat)
    pheno_tbl[["progression"]] <- sample(
        c(0, 1),
        size = n_samples,
        replace = TRUE
    )
    pheno_tbl[["pfs_years"]] <- rnorm(n_samples, 2, 1)
    pheno_tbl[["discrete_var"]] <- sample(1:3, size = n_samples, replace = TRUE)
    pheno_tbl[["continuous_var"]] <- rnorm(n_samples, 10, 10)
    pheno_tbl[["ipi"]] <- sample(1:5, size = n_samples, replace = TRUE)
    # insert NAs
    na_rows <- sample(1:n_samples, n_na_in_pheno, replace = TRUE)
    na_cols <- sample(2:ncol(pheno_tbl), n_na_in_pheno, replace = TRUE)
    if(n_na_in_pheno > 0){
        for(i in 1:n_na_in_pheno){
            pheno_tbl[na_rows[i], na_cols[i]] <- NA
        }
    }

    if(is.character(to_csv)){
        readr::write_csv(expr_tbl, file.path(to_csv, "expr.csv"))
        readr::write_csv(pheno_tbl, file.path(to_csv, "pheno.csv"))
    }

    res <- list(
        "expr_mat" = expr_mat,
        "pheno_tbl" = pheno_tbl
    )
    return(res)
}