new_DataSpec <- function(
    patient_id_col,
    pfs_col,
    progression_col,
    gene_id_col
){
    stopifnot(is.character(patient_id_col))
    stopifnot(is.character(pfs_col))
    stopifnot(is.character(progression_col))
    stopifnot(is.character(gene_id_col))
    data_spec_list <- list(
        "patient_id_col" = patient_id_col,
        "pfs_col" = pfs_col,
        "progression_col" = progression_col,
        "gene_id_col" = gene_id_col
    )
    return(structure(data_spec_list, class = "DataSpec")) 
}

#' @title Construct a DataSpec S3 object
#' @description A DataSpec object holds the pieces of information key to reading and
#' preparing the expression and pheno data, namely:
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
#' @details Note that the appearance of the pheno data changes over the course of the
#' pipeline: first, it is csv file, later it is a tibble.
#' @export
DataSpec <- function(
    patient_id_col = "patient_id",
    pfs_col = "pfs_years",
    progression_col = "progression",
    gene_id_col = "gene_id"
){
    data_spec <- new_DataSpec(
        patient_id_col = patient_id_col,
        pfs_col = pfs_col,
        progression_col = progression_col,
        gene_id_col = gene_id_col
    )
    return(data_spec)
}