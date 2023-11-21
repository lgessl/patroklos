#' @title Read expression and pheno data from csv files into consistent format
#' @description Read expression data into a matrix and pheno data into a tibble.
#' Both will hold patients as rows.
#' @param directory string. The directory where both expression and pheno csv
#' files lie.
#' @param expr_fname string. The name of the expression .csv file inside `directory`.
#' See details for the expected format.
#' @param pheno_fname string. The name of the pheno data .csv inside `directory`.
#' See details for the expected format.
#' @param patient_id_col string. The name of the column in the pheno data file
#' that holds the patient identifiers.
#' @param gene_id_col string. The name of the column in the expression data file
#' that holds the gene identifiers.
#' @return A list with a numeric matrix, named `expr`, and a tibble named `pheno`. 
#' `expr` holds the expression data, with patient ids as row names and gene ids as
#' column names. I.e., we transpose the expression data. `pheno` holds the pheno data, 
#' with the patient ids in the first column `patient_id_col`.
#' @details The pheno csv file holds the samples as rows (with the unique sample names
#'  in the column `patient_id`), the variables as columns. We need at least the columns
#'  * progression (integer, 0 for censored, 1 for uncensored),
#'  * pfs_yrs (integer, the time in years to progression or censoring),
#'  * `ìnclude_from_continuous_pheno`,
#'  * `include_from_discrete_pheno` and
#'  * `patient_id_col`.
#' The expr csv file holds the genes as rows (with the gene ids in a row called 
#' `gene_id_col`), the samples as columns. The sample names in pheno and expression files
#' must be identical and in the same order.
read <- function(
    directory,
    expr_fname,
    pheno_fname,
    patient_id_col,
    gene_id_col
){
    # read
    fnames <- c(expr_fname, pheno_fname)
    res <- list("expr" = NULL, "pheno" = NULL)
    for (i in 1:length(fnames)){
        full_path <- file.path(directory, fnames[i])
        res[[i]] <- readr::read_csv(full_path, show_col_types = FALSE)
    }
    expr_tbl <- res[["expr"]]
    pheno_tbl <- res[["pheno"]]

    # expression to matrix
    gene_names <- res[["expr"]][[gene_id_col]]
    expr_mat <- expr_tbl |>
        dplyr::select(!all_of(gene_id_col)) |>
        as.matrix() |> 
        t()
    colnames(expr_mat) <- gene_names

    # pheno: move patient ids into first column
    pheno_tbl <- pheno_tbl |>
        dplyr::relocate(all_of(patient_id_col))

    return(res)
}