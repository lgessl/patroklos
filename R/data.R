#' @title An R6 class for the data
#' @description A `Data` object holds the phenotype and expression data belonging to
#' a data set. It specifies the names of important features in the columns, reads in the data, puts 
#' the spotlight on a part of the data: the cohort, and prepares the data for a model.
#' @export
Data <- R6::R6Class("Data",
    public = list(

        #' @field name A telling name for the data set.
        name = NULL,
        #' @field directory Directory where the expression and pheno csv files lie.
        directory = NULL,
        #' @field pivot_time_cutoff Time cutoff that divides the sample into low-risk (event 
        #' before) and high-risk (event after).
        #' assessment.
        pivot_time_cutoff = NULL,
        #' @field cohort Regular expression to subset the data to a cohort.
        cohort = NULL,
        #' @field imputer Function handling NAs in the predictor matrix.
        imputer = NULL,
        #' @field expr_mat Named numeric matrix. Samples correspond to rows.
        expr_mat = NULL,
        #' @field pheno_tbl A tibble with phenotypic features and samples as rows.
        pheno_tbl = NULL,
        #' @field expr_file Name of the expression csv file inside `directory`.
        expr_file = NULL,
        #' @field pheno_file Name of the pheno data csv inside `directory`.
        pheno_file = NULL,
        #' @field cohort_col Find the cohort of a sample in this column of the
        #' pheno data.
        cohort_col = NULL,
        #' @field patient_id_col The name of the column in the pheno data that holds unique
        #' patient identifiers.
        patient_id_col = NULL,
        #' @field time_to_event_col The name of the column in the pheno data that holds the
        #' time-to-event values.
        time_to_event_col = NULL,
        #' @field event_col The name of the column in the pheno data that holds
        #' the event status encoded as 1 = occurrence, 0 = censoring.
        event_col = NULL,
        #' @field gene_id_col The name of the column in the expression data that holds
        #' the gene identifiers.
        gene_id_col = NULL,
        #' @field benchmark_col The name of the column in the pheno data that holds the
        #' benchmark risk score (like the IPI).
        benchmark_col = NULL,

        #' @description Construct a `Data` R6 object.
        #' @param name string. A telling name for the data set.
        #' @param directory string. The directory where both expression and pheno csv
        #' files lie.
        #' @param pivot_time_cutoff numeric. Time cutoff that divides the samples into low-risk 
        #' (event before) and high-risk (event after).
        #' @param cohort string. At the end of preparing the data, subset it to those
        #' samples whose value in the `cohort_col` column matches `cohort`.
        #' @param imputer function or `NULL`. Function imputing `NA`s in the predictor matrix.
        #' See [`imputer_prototype()`] for its interface. Default is [`mean_impute()`]. `NULL` 
        #' means no imputation.
        #' @param time_to_event_col string. The name of the column in the pheno data that holds the
        #' time-to-event values.
        #' @param event_col string. The name of the column in the pheno data that holds
        #' the event status encoded as 1 = occurrence, 0 = censoring.
        #' @param cohort_col string. The name of the column in the pheno data that holds
        #' the cohort a sample belongs to.
        #' @param benchmark_col string or `NULL`. The name of the column in the pheno data that
        #' holds the output of a benchmark model.
        #' @param expr_file string. The name of the expression csv file inside `directory`.
        #' Default is `"expr.csv"`. See details for the expected format.
        #' @param pheno_file string. The name of the pheno data csv inside `directory`.
        #' Default is `"pheno.csv"`. See details for the expected format.
        #' @param patient_id_col string. The name of the column in the pheno data that holds
        #' the patient identifiers.
        #' @param gene_id_col string. The name of the column in the expression data that holds
        #' the gene identifiers.
        #' @return A Data object.
        #' @details The pheno csv file holds the samples as rows (with *unique* sample ids in the
        #' first (character) column called `patient_id_col`), the variables as columns.
        #'
        #' The expression csv file holds the genes as rows (with *unique* gene ids in the first
        #' (character) column called `gene_id_col`), the samples as columns.
        initialize = function(
            name,
            directory,
            pivot_time_cutoff,
            cohort,
            imputer = mean_impute,
            time_to_event_col,
            event_col,
            cohort_col,
            benchmark_col = NULL,
            expr_file = "expr.csv",
            pheno_file = "pheno.csv",
            patient_id_col = "patient_id",
            gene_id_col = "gene_id"
        ) {
            data_initialize(
                self, private, name, directory, pivot_time_cutoff,
                expr_file, pheno_file, cohort, patient_id_col, time_to_event_col,
                event_col, cohort_col, benchmark_col, gene_id_col, imputer
            )
        },

        #' @description Read expression data into the `expr_mat` attribute and pheno data into the 
        #' `pheno_tbl` attribute.
        read = function() {
            data_read(self, private)
        },

        #' @description Prepare the already read-in data for a model.
        #' @param model A `Model` object.
        #' @param quiet logical. If TRUE, suppress messages.
        #' @details You need to set `cohort` before calling this method.
        prepare = function(model, quiet = FALSE) {
            data_prepare(self, private, model, quiet)
        },

        #' @description Calculate the quantiles of the survival times.
        #' @param round_digits integer. Round the numbers in the returned tibble
        #' to this number of digits after the point.
        #' @return A tibble with two columns. For each quantile q in the first
        #' column, the time-to-event value in the second column.
        #' @details We take censoring into account.
        survival_quantiles = function(round_digits = 3) {
            data_survival_quantiles(self, private, round_digits)
        },

        #' @description Split the data into a train and test cohort
        #' @param train_prop numeric. Proportion of the data to put in the train 
        #' cohort.
        #' @param save logical. If TRUE, save the named cohort vector to a file.
        #' @param keep_risk logical. If TRUE, keep the ratio of high-risk versus 
        #' low-risk samples in train and test cohort the same as in the complete
        #' data set. 
        #' @param quiet logical. If TRUE, suppress messages.
        #' @details Cohort affiliation will show up in the column `cohort_col` in `pheno_tbl`
        split = function(train_prop, save = TRUE, keep_risk = TRUE, quiet = FALSE)
            data_split(self, private, train_prop, save, keep_risk, quiet),

        #' @description Quality control at the end of preprocessing
        #' @param expr_tbl A tibble with the expression data, the first column, named `gene_id_col`,
        #' holds the gene identifiers and the other columns the samples.
        #' @details Check if the expression and pheno tibble are consistent
        #' with the other attributes of the Data object. *You* typically call this
        #' method at the end of preprocessing, and the read() method calls it.
        qc_preprocess = function(expr_tbl) {
            data_qc_preprocess(self, private, expr_tbl)
        }
    )
)

#' @import R6
NULL