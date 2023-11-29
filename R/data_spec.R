new_DataSpec <- function(
    directory,
    expr_fname,
    pheno_fname,
    patient_id_col,
    pfs_col,
    progression_col,
    gene_id_col
){
    stopifnot(is.character(directory))
    stopifnot(is.character(expr_fname))
    stopifnot(is.character(pheno_fname))
    stopifnot(is.character(patient_id_col))
    stopifnot(is.character(pfs_col))
    stopifnot(is.character(progression_col))
    stopifnot(is.character(gene_id_col))
    data_spec_list <- list(
        "directory" = directory,
        "expr_fname" = expr_fname,
        "pheno_fname" = pheno_fname,
        "patient_id_col" = patient_id_col,
        "pfs_col" = pfs_col,
        "progression_col" = progression_col,
        "gene_id_col" = gene_id_col
    )
    return(structure(data_spec_list, class = "DataSpec")) 
}

#' @title Construct a DataSpec S3 object
#' @description A DataSpec object specifies the location and format of the expression
#' and pheno data of a single data set. This enables reading and preparing the data.
#' @param directory string. The directory where both expression and pheno csv
#' files lie.
#' @param expr_fname string. The name of the expression csv file inside `directory`.
#' Default is `"expr.csv"`. See details for the expected format.
#' @param pheno_fname string. The name of the pheno data csv inside `directory`.
#' Default is `"pheno.csv"`. See details for the expected format.
#' @param patient_id_col string. The name of the column in the pheno data that holds 
#' the patient identifiers. Default is `"patient_id"`.
#' @param pfs_col string. The name of the column in the pheno data that holds the 
#' progression-free survival (PFS) values. Default is `"pfs_years"`.
#' @param progression_col string. The name of the column in the pheno data that holds
#' the progression status encoded as 1 = progression, 0 = no progression. Default is
#' `"progression"`.
#' @param gene_id_col string. The name of the column in the expression data that holds
#' the gene identifiers. Default is `"gene_id"`.
#' @return A DataSpec object.
#' @details The pheno csv file holds the samples as rows (with *unique* sample ids in a
#' column called `patient_id_col`), the variables as columns. The expr csv file holds 
#' the genes as rows (with *unique* gene ids in a column called `gene_id_col`), the 
#' samples as columns.
#' While the computatinal representation of both expression and pheno data will
#' change over the course of the pipeline, a DataSpec object will hold timeless
#' information on the data.
#' @export
DataSpec <- function(
    directory = ".",
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    patient_id_col = "patient_id",
    pfs_col = "pfs_years",
    progression_col = "progression",
    gene_id_col = "gene_id"
){
    data_spec <- new_DataSpec(
        directory = directory,
        expr_fname = expr_fname,
        pheno_fname = pheno_fname,
        patient_id_col = patient_id_col,
        pfs_col = pfs_col,
        progression_col = progression_col,
        gene_id_col = gene_id_col
    )
    return(data_spec)
}