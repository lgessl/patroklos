#' @title Subset existing data to mock data
#' @description Read in the expression and pheno data of an existing dataset,
#' randomly pick a subset of samples and genes, and return and possibly save
#' the subset as mock data.
#' @param directory string. The directory where both expression and pheno csv
#' files lie.
#' @param n_samples integer. The number of samples to randomly pick from the
#' initial data.
#' @param n_genes integer. The number of genes to randomly pick from the
#' initial data.
#' @param expr_fname string. The name of the expression csv file inside `directory`.
#' See details of `read` for the expected format. Default is `"expr.csv"`.
#' @param pheno_fname string. The name of the pheno data csv inside `directory`.
#' See details of `read` for the expected format. Default is `"pheno.csv"`.
#' @param patient_id_col string. The name of the column in the pheno data csv
#' file that holds the patient identifiers. Default is `"patient_id"`.
#' @param gene_id_col string. The name of the column in the expression data csv
#' file that holds the gene identifiers. Default is `"gene_id"`.
#' @param save logical. Whether to save the subset as csv files. Default is 
#' `TRUE`.
#' @param save_suffix string. The suffix to append to the file names of the
#' subset csv files. Default is `"mock"`.
#' @return A list with two tibbles, `expr` and `pheno`, the subset expression
#' and pheno data. If `save` is `TRUE`, the subset data will be saved as csv
#' files in the format desired by `read`.
mock_data_from_existing <- function(
    directory,
    n_samples = 50,
    n_genes = 10,
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    patient_id_col = "patient_id",
    gene_id_col = "gene_id",
    save = TRUE,
    save_suffix = "mock"
){
    expr_tbl <- readr::read_csv(
        file.path(directory, expr_fname),
        show_col_types = FALSE
    )
    pheno_tbl <- readr::read_csv(
        file.path(directory, pheno_fname),
        show_col_types = FALSE
    )
    # randomly pick
    patients_sub <- sample(
        pheno_tbl[[patient_id_col]],
        size = n_samples,
        replace = FALSE
        )
    genes_sub <- sample(
        expr_tbl[[gene_id_col]],
        size = n_genes,
        replace = FALSE
        )
    # subset
    expr_tbl <- expr_tbl[
        expr_tbl[[gene_id_col]] %in% genes_sub,
        c(gene_id_col, patients_sub) # first column are gene names
        ]
    pheno_tbl <- pheno_tbl[pheno_tbl[[patient_id_col]] %in% patients_sub, ]
    # save
    if(save){
        expr_fname <- stringr::str_replace(
            expr_fname, 
            ".csv", 
            stringr::str_c("_", save_suffix, ".csv")
        )
        pheno_fname <- stringr::str_replace(
            pheno_fname, 
            ".csv", 
            stringr::str_c("_", save_suffix, ".csv")
        )
        readr::write_csv(
            expr_tbl, 
            file.path(directory, expr_fname)
        )
        readr::write_csv(
            pheno_tbl, 
            file.path(directory, pheno_fname)
        )
    }
    res <- list("expr" = expr_tbl, "pheno" = pheno_tbl)
    return(res)
}