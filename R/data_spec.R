new_DataSpec <- function(
    name,
    directory,
    expr_fname,
    pheno_fname,
    patient_id_col,
    pfs_col,
    progression_col,
    ipi_col,
    ipi_feat_cols,
    gene_id_col
){
    stopifnot(is.character(name))
    stopifnot(is.character(directory))
    stopifnot(is.character(expr_fname))
    stopifnot(is.character(pheno_fname))
    stopifnot(is.character(patient_id_col))
    stopifnot(is.character(pfs_col))
    stopifnot(is.character(progression_col))
    stopifnot(is.character(ipi_col))
    stopifnot(is.character(ipi_feat_cols) | is.null(ipi_feat_cols))
    stopifnot(is.character(gene_id_col))

    # More complex checks for ipi_feat_cols
    if(is.character(ipi_feat_cols)){
        if(is.null(names(ipi_feat_cols))){
            stop("ipi_feat_cols must be a named character vector.")
        }
        ipi_feat_names <- names(ipi_feat_cols_default)
        if(any(names(ipi_feat_cols) != ipi_feat_names)){
            stop("ipi_feat_cols must be a character vector with the following names: ",
                paste0(ipi_feat_names, collapse = ", "))
        }
    }

    data_spec_list <- list(
        "name" = name,
        "directory" = directory,
        "expr_fname" = expr_fname,
        "pheno_fname" = pheno_fname,
        "patient_id_col" = patient_id_col,
        "pfs_col" = pfs_col,
        "progression_col" = progression_col,
        "ipi_col" = ipi_col,
        "ipi_feat_cols" = ipi_feat_cols,
        "gene_id_col" = gene_id_col
    )
    return(structure(data_spec_list, class = "DataSpec")) 
}


# For below helper 
ipi_feat_cols_default <- c(
    "age" = "age",
    "stage" = "ann_arbor_stage",
    "ldh_ratio" = "ldh_ratio",
    "performance_status" = "ecog_performance_status",
    "n_extranodal_sites" = "n_extranodal_sites"
)


#' @title Construct a DataSpec S3 object
#' @description A DataSpec object specifies the location and format of the expression
#' and pheno data of a single data set. This enables reading and preparing the data.
#' @param name string. A telling name for the data set (e.g. `"Schmitz et al. 2018"`).
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
#' @param ipi_col string. The name of the column in the pheno data that holds the
#' International Prognostic Index (IPI) values. Default is `"ipi"`.
#' @param ipi_feat_cols named character vector of length 5 or `NULL`. The names of the 
#' *five* features in the pheno data needed to compute the IPI. If NULL, no information 
#' about the IPI features is provided, e.g. because they are not in the pheno data. For 
#' the default, see Usage.
#' @param gene_id_col string. The name of the column in the expression data that holds
#' the gene identifiers. Default is `"gene_id"`.
#' @return A DataSpec object.
#' @details The pheno csv file holds the samples as rows (with *unique* sample ids in the
#' first column called `patient_id_col`), the variables as columns. The expr csv file 
#' holds the genes as rows (with *unique* gene ids in the first column called `gene_id_col`), 
#' the samples as columns.
#' While the computational representation of both expression and pheno data will
#' change over the course of the pipeline, a DataSpec object will hold timeless
#' information on the data.
#' @export
DataSpec <- function(
    name,
    directory = ".",
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    patient_id_col = "patient_id",
    pfs_col = "pfs_years",
    progression_col = "progression",
    ipi_col = "ipi",
    ipi_feat_cols = c(
        "age" = "age",
        "stage" = "ann_arbor_stage",
        "ldh_ratio" = "ldh_ratio",
        "performance_status" = "ecog_performance_status",
        "n_extranodal_sites" = "n_extranodal_sites"
    ),
    gene_id_col = "gene_id"
){
    data_spec <- new_DataSpec(
        name = name,
        directory = directory,
        expr_fname = expr_fname,
        pheno_fname = pheno_fname,
        patient_id_col = patient_id_col,
        pfs_col = pfs_col,
        progression_col = progression_col,
        ipi_col = ipi_col,
        ipi_feat_cols = ipi_feat_cols,
        gene_id_col = gene_id_col
    )
    return(data_spec)
}


#' @title To every `ModelSpec` in a list of `ModelSpec`s, prepend a fixed directory to the
#' `save_dir` attribute
#' @description Often one specifies the models in general, for all data sets. If you fit the 
#' models to a specific data set (say `"mock"`), you might want to prepend a fixed directory 
#' like `"results/mock"` to the `save_dir` attribute of all `ModelSpec`s in the list.
#' @param model_spec_list list of `ModelSpec`s.
#' @param dir string. The directory to prepend to the `save_dir` attribute of all 
#' `ModelSpec`s in `model_spec_list`.
#' @return A list of `ModelSpec`s, with `dir` prepended to the `save_dir` attribute.
#' @export
prepend_to_save_dir <- function(
    model_spec_list, 
    dir
){
    stopifnot(is.character(dir))
    for(i in seq_along(model_spec_list)){
        model_spec_list[[i]]$save_dir <- file.path(dir, model_spec_list[[i]]$save_dir)
    }
    return(model_spec_list)
}