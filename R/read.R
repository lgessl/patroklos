#' @title Read expression and pheno data from csv files into tibbles
#' @description Read expression and pheno csv files into tibbles and move
#' patient and gene identifiers to row names.
#' @param directory string. The directory where both expression and pheno csv
#' files lie.
#' @param expr_fname string. The name of the expression .csv file inside `directory`.
#' See details for the expected format.
#' @param pheno_fname string. The name of the pheno data .csv inside `directory`.
#' See details for the expected format.
#' @param patient_id_col string. The name of the column in the pheno data file
#' that holds the patient identifiers.
#' @param gene_col string. The name of the column in the expression data file
#' that holds the gene identifiers.
#' @details The pheno .csv files holds the samples as rows (with the unique sample names
#'  in the very first column), the variables as columns. We need at least the columns
#'  * progression (integer, 0 for censored, 1 for uncensored),
#'  * pfs_yrs (integer, the time in years to progression or censoring),
#'  * `ìnclude_from_continuous_pheno`,
#'  * `include_from_discrete_pheno` and
#'  * `patient_id_col`.
#' The expr .csv file holds the genes as rows (with the unique genes in a row called 
#' `gene_col`), the samples as columns. The sample names in pheno and expression files
#' must be identical and in the same order.
read <- function(
    directory,
    expr_fname,
    pheno_fname,
    patient_id_col,
    gene_col
){
    fnames <- c(expr_fname, pheno_fname)
    res <- list("expr" = NULL, "pheno" = NULL)
    for (i in 1:length(fnames)){
        full_path <- file.path(data_dir, fnames[i])
        res[[i]] <- readr::read_csv(full_path)
    }
    
    # actual data is actual data frame
    res[["expr"]] <- expr_df %>% 
        tibble::column_to_rownames(var = "Gene")
    res[["pheno"]] <- pheno_df %>%
        tibble::column_to_rownames(var = "patient_id")

    return(res)
}