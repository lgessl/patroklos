#' @title An R6 class for the data
#' @description A Data object holds the phenotype and expression data belonging to 
#' a data set. It specifies the names of important features in the columns, reads from 
#' and writes to csv files and takes care of splitting the data into train and test 
#' cohort.
Data <- R6::R6Class("Data",
    public = list(

        #' @field name How to refer to the data in plots and messages.
        name = NULL,
        #' @field directory Directory where the expression and pheno csv files lie.
        directory = NULL,
        #' @field train_prop Proportion of samples in the training cohort.
        train_prop = NULL,
        #' @field expr_mat The expression matrix.
        expr_mat = NULL,
        #' @field pheno_tbl The phenotype tibble.
        pheno_tbl = NULL,
        #' @field expr_file Name of the expression csv file inside `directory`.
        expr_file = NULL,
        #' @field pheno_file Name of the pheno data csv inside `directory`.
        pheno_file = NULL,
        #' @field cohort The cohort of the data, train or test.
        cohort = NULL,
        #' @field patient_id_col The name of the column in the pheno data that holds
        #' the patient identifiers.
        patient_id_col = NULL,
        #' @field time_to_event_col The name of the column in the pheno data that holds the
        #' time-to-event values.
        time_to_event_col = NULL,
        #' @field event_col The name of the column in the pheno data that holds
        #' the event status encoded as 1 = occurrence, 0 = censoring.
        event_col = NULL,
        #' @field benchmark_col The name of the column in the pheno data that holds the
        #' benchmark risk score (like the IPI).
        benchmark_col = NULL,
        #' @field gene_id_col The name of the column in the expression data that holds
        #' the gene identifiers.
        gene_id_col = NULL,
        #' @field split_col_prefix Column-name prefix for those columns holding splits into
        #' train and test cohort.
        split_col_prefix = NULL,
        #' @field pivot_time_cutoff Pivotal time cutoff for the analysis and, explicitly,
        #' for splitting the data into the train and test cohort.
        pivot_time_cutoff = NULL,

        #' @description A Data object specifies the location and format of the expression
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
        #' @param expr_file string. The name of the expression csv file inside `directory`.
        #' Default is `"expr.csv"`. See details for the expected format.
        #' @param pheno_file string. The name of the pheno data csv inside `directory`.
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
        #' @return A Data object.
        #' @details The pheno csv file holds the samples as rows (with *unique* sample ids in the
        #' first column called `patient_id_col`), the variables as columns. The expr csv file 
        #' holds the genes as rows (with *unique* gene ids in the first column called `gene_id_col`), 
        #' the samples as columns.
        #' While the computational representation of both expression and pheno data will
        #' change over the course of the pipeline, a Data object will hold timeless
        #' information on the data.
        initialize = function(
            name,
            directory,
            train_prop,
            pivot_time_cutoff = 2.0,
            expr_file = "expr.csv",
            pheno_file = "pheno.csv",
            cohort = NULL,
            patient_id_col = "patient_id",
            time_to_event_col = "pfs_years",
            event_col = "progression",
            split_col_prefix = "split_",
            benchmark_col = NULL,
            gene_id_col = "gene_id"
        )
            data_initialize(self, private, name, directory, train_prop,              
                pivot_time_cutoff, expr_file, pheno_file, cohort, patient_id_col, 
                time_to_event_col, split_col_prefix, event_col, benchmark_col, 
                gene_id_col),

        #' @description Read expression data into the matrix `expr_mat` and pheno 
        #' data into the tibble `pheno_tbl`. Both will hold patients as rows.
        read = function()
            data_read(self, private),

        #' @description Prepare the read-in data for a model.
        #' @param model A Model object.
        #' @details You need to set the cohort before calling this method.
        prepare = function(model)
            data_prepare(self, private, model)
    ),

    private = list(
        dummy = "dummy"
    )
)

#' @import R6
NULL