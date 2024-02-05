new_DataSpec <- function(
    name,
    directory,
    train_prop,
    pivot_time_cutoff,
    expr_fname,
    pheno_fname,
    cohort,
    patient_id_col,
    time_to_event_col,
    event_col,
    benchmark_col,
    gene_id_col,
    split_col_prefix
){
    stopifnot(is.character(name))
    stopifnot(is.character(directory))
    stopifnot(is.numeric(train_prop))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(is.character(expr_fname))
    stopifnot(is.character(pheno_fname))
    stopifnot(is.null(cohort) || is.character(cohort))
    stopifnot(is.character(patient_id_col))
    stopifnot(is.character(time_to_event_col))
    stopifnot(is.character(event_col))
    stopifnot(is.character(benchmark_col) || is.null(benchmark_col))
    stopifnot(is.character(gene_id_col))
    stopifnot(is.character(split_col_prefix))

    data_spec_list <- list(
        "name" = name,
        "directory" = directory,
        "train_prop" = train_prop,
        "pivot_time_cutoff" = pivot_time_cutoff,
        "expr_fname" = expr_fname,
        "pheno_fname" = pheno_fname,
        "cohort" = cohort,
        "patient_id_col" = patient_id_col,
        "time_to_event_col" = time_to_event_col,
        "event_col" = event_col,
        "benchmark_col" = benchmark_col,
        "gene_id_col" = gene_id_col,
        "split_col_prefix" = split_col_prefix
    )
    return(structure(data_spec_list, class = "DataSpec")) 
}


#' @title Construct a DataSpec S3 object
#' @description A DataSpec object specifies the location and format of the expression
#' and pheno data of a single data set. This enables reading and preparing the data.
#' @param name string. A telling name for the data set (e.g. `"Schmitz et al. 2018"`).
#' @param directory string. The directory where both expression and pheno csv
#' files lie.
#' @param train_prop numeric in (0, 1). The proportion of samples to draw for the 
#' training cohort in every splitting run.
#' @param pivot_time_cutoff numeric or NULL. Pivotal time cutoff for the analysis and, explicitly,
#' for splitting the data into the train and test cohort: If a numeric in (0, 1), preserve 
#' the the proportion of individuals below and above `time_cutoff` in both cohorts. Default 
#' is `NULL`, i.e. no further constraints on splitting.
#' @param expr_fname string. The name of the expression csv file inside `directory`.
#' Default is `"expr.csv"`. See details for the expected format.
#' @param pheno_fname string. The name of the pheno data csv inside `directory`.
#' Default is `"pheno.csv"`. See details for the expected format.
#' @param cohort string in `c("train", "test")` or `NULL`. The cohort, train or test, to 
#' prepare the data for. If NULL, the default, some functions will set it themselves
#' (e.g. `training_camp()` to `"train"`, `assess_2d()` to `"test"`), others will
#' throw an error.
#' @param patient_id_col string. The name of the column in the pheno data that holds 
#' the patient identifiers. Default is `"patient_id"`.
#' @param time_to_event_col string. The name of the column in the pheno data that holds the 
#' time-to-event values. Default is `"pfs_years"`.
#' @param event_col string. The name of the column in the pheno data that holds
#' the event status encoded as 1 = occurrence, 0 = censoring. Default is
#' `"progression"`.
#' @param split_col_prefix string. Column-name prefix for those columns holding splits into 
#' train and test cohort. Some of these columns may be present in the pheno data already,
#' others will be added during the training process if required. Default is `"split_"`, 
#' i.e., split columns are named `"split_1"`, `"split_2"`, etc.
#' @param benchmark_col string or `NULL`. The name of the column in the pheno data that 
#' holds the benchmark risk score (like the IPI). Default is `NULL`.
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
    directory,
    train_prop,
    pivot_time_cutoff = 2.0,
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    cohort = NULL,
    patient_id_col = "patient_id",
    time_to_event_col = "pfs_years",
    split_col_prefix = "split_",
    event_col = "progression",
    benchmark_col = NULL,
    gene_id_col = "gene_id"
){
    if(!is.null(cohort)){
        cohort <- match.arg(cohort, c("train", "test"))
    }
    data_spec <- new_DataSpec(
        name = name,
        directory = directory,
        train_prop = train_prop,
        pivot_time_cutoff = pivot_time_cutoff,
        expr_fname = expr_fname,
        pheno_fname = pheno_fname,
        cohort = cohort,
        patient_id_col = patient_id_col,
        time_to_event_col = time_to_event_col,
        event_col = event_col,
        benchmark_col = benchmark_col,
        gene_id_col = gene_id_col,
        split_col_prefix = split_col_prefix
    )
    return(data_spec)
}


#' @title To every `ModelSpec` in a list of `ModelSpec`s, prepend a fixed directory to the
#' `directory` attribute
#' @description Often one specifies the models in general, for all data sets. If you fit the 
#' models to a specific data set (say `"mock"`), you might want to prepend a fixed directory 
#' like `"results/mock"` to the `directory` attribute of all `ModelSpec`s in the list.
#' @param model_spec_list list of `ModelSpec`s.
#' @param dir string. The directory to prepend to the `directory` attribute of all 
#' `ModelSpec`s in `model_spec_list`.
#' @return A list of `ModelSpec`s, with `dir` prepended to the `directory` attribute.
#' @export
prepend_to_directory <- function(
    model_spec_list, 
    dir
){
    stopifnot(is.character(dir))
    for(i in seq_along(model_spec_list)){
        model_spec_list[[i]]$directory <- file.path(dir, model_spec_list[[i]]$directory)
    }
    return(model_spec_list)
}