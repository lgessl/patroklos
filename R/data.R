#' @title An R6 class for the data
#' @description A Data object holds the phenotype and expression data belonging to 
#' a data set. It specifies the names of important features in the columns, reads from 
#' and writes to csv files and puts the spotlight on a part of the data: the 
#' cohort.
#' @export
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
        #' @field cohort_col Find the cohort of a sample in this column of the 
        #' pheno data.
        cohort_col = NULL,
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
        #' @field pivot_time_cutoff Pivotal time cutoff for the analysis and 
        #' assessment.
        pivot_time_cutoff = NULL,
        #' @field imputer Function handling NAs in the predictor matrix.
        imputer = NULL,

        #' @description A Data object specifies the location and format of the expression
        #' and pheno data of a single data set. This enables reading and preparing the data.
        #' @param name string. A telling name for the data set (e.g. `"Schmitz et al. 2018"`).
        #' @param directory string. The directory where both expression and pheno csv
        #' files lie.
        #' @param pivot_time_cutoff numeric. Pivotal time cutoff for the analysis.
        #' @param cohort string. At the end of preparing the data, subset it to those 
        #' samples whose value in the `cohort_col` column matches `cohort`. 
        #' @param time_to_event_col string. The name of the column in the pheno data that holds the 
        #' time-to-event values.
        #' @param event_col string. The name of the column in the pheno data that holds
        #' the event status encoded as 1 = occurrence, 0 = censoring.
        #' @param cohort_col string. The name of the column in the pheno data that holds
        #' the cohort a sample belongs to.
        #' @param benchmark_col string or `NULL`. The name of the column in the pheno data that 
        #' holds the benchmark risk score (like the IPI).
        #' @param expr_file string. The name of the expression csv file inside `directory`.
        #' Default is `"expr.csv"`. See details for the expected format.
        #' @param pheno_file string. The name of the pheno data csv inside `directory`.
        #' Default is `"pheno.csv"`. See details for the expected format.
        #' @param patient_id_col string. The name of the column in the pheno data that holds 
        #' the patient identifiers.
        #' @param gene_id_col string. The name of the column in the expression data that holds
        #' the gene identifiers.
        #' @param imputer function or NULL. Function handling NAs in the predictor matrix.
        #' It takes a single argument `x`, the predictor matrix, a numeric with categorical 
        #' variables dichotomized to binary dummy variables, and returns a matrix of the 
        #' same shape with non-NAs left untouched. Default is [`mean_impute()`]. 
        #' NULL means no imputation (or applying the identity function on `x`).
        #' @return A Data object.
        #' @details The pheno csv file holds the samples as rows (with *unique* sample ids in the
        #' first (character) column called `patient_id_col`), the variables as columns. 
        #' 
        #' The expression csv file holds the genes as rows (with *unique* gene ids in the first 
        #' (character) column called `gene_id_col`), the samples as columns.
        #' 
        #' While the computational representation of both expression and pheno data will
        #' change over the course of the pipeline, a Data object will hold timeless
        #' information on the data.
        initialize = function(
            name,
            directory,
            pivot_time_cutoff,
            time_to_event_col,
            cohort,
            event_col = "progression",
            cohort_col = "split_1",
            benchmark_col = NULL,
            expr_file = "expr.csv",
            pheno_file = "pheno.csv",
            patient_id_col = "patient_id",
            gene_id_col = "gene_id",
            imputer = mean_impute
        )
            data_initialize(self, private, name, directory, train_prop,              
                pivot_time_cutoff, expr_file, pheno_file, cohort, patient_id_col, 
                time_to_event_col, event_col, cohort_col, benchmark_col, 
                gene_id_col, imputer),

        #' @description Read expression data into the matrix `expr_mat` and pheno 
        #' data into the tibble `pheno_tbl`. Both will hold patients as rows.
        read = function()
            data_read(self, private),

        #' @description Prepare the read-in data for a model.
        #' @param model A Model object.
        #' @param quiet logical. If TRUE, suppress messages.
        #' @details You need to set the cohort before calling this method.
        prepare = function(model, quiet = FALSE)
            data_prepare(self, private, model, quiet),

        #' @description Calculate the quantiles of the survival times.
        #' @param round_digits integer. Round the numbers in the returned tibble 
        #' to this number of digits after the point.
        #' @return A tibble with two columns. For each quantile q in the first 
        #' column, the time-to-event value in the second column.
        #' @details We take censoring into account. 
        survival_quantiles = function(round_digits = 3)
            data_survival_quantiles(self, private, round_digits),

        #' @description Quality control at the end of preprocessing
        #' @param expr_tbl A tibble with the expression data, the product of 
        #' reading in an expression csv file as described in the details of 
        #' the initialize() method.
        #' @details Check if the expression and pheno tibble are consistent 
        #' with the other attributes of the Data object. You typically call this 
        #' method at the end of preprocessing, and the read() method calls it. 
        qc_preprocess = function(expr_tbl)
            data_qc_preprocess(self, private, expr_tbl)
    )
)

#' @import R6
NULL