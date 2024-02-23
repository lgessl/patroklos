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
        #' @return A DataSpec object.
        #' @details The pheno csv file holds the samples as rows (with *unique* sample ids in the
        #' first column called `patient_id_col`), the variables as columns. The expr csv file 
        #' holds the genes as rows (with *unique* gene ids in the first column called `gene_id_col`), 
        #' the samples as columns.
        #' While the computational representation of both expression and pheno data will
        #' change over the course of the pipeline, a DataSpec object will hold timeless
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
            split_col_prefix = "split_",
            event_col = "progression",
            benchmark_col = NULL,
            gene_id_col = "gene_id"
        ){
            if(!is.null(cohort)){
                cohort <- match.arg(cohort, c("train", "test"))
            }
            stopifnot(is.character(name))
            stopifnot(is.character(directory))
            stopifnot(is.numeric(train_prop))
            stopifnot(is.numeric(pivot_time_cutoff))
            stopifnot(is.character(expr_file))
            stopifnot(is.character(pheno_file))
            stopifnot(is.null(cohort) || is.character(cohort))
            stopifnot(is.character(patient_id_col))
            stopifnot(is.character(time_to_event_col))
            stopifnot(is.character(event_col))
            stopifnot(is.character(benchmark_col) || is.null(benchmark_col))
            stopifnot(is.character(gene_id_col))
            stopifnot(is.character(split_col_prefix))

            self$name <- name
            self$directory <- directory
            self$train_prop <- train_prop
            self$pivot_time_cutoff <- pivot_time_cutoff
            self$expr_file <- expr_file
            self$pheno_file <- pheno_file
            self$cohort <- cohort
            self$patient_id_col <- patient_id_col
            self$time_to_event_col <- time_to_event_col
            self$event_col <- event_col
            self$benchmark_col <- benchmark_col
            self$gene_id_col <- gene_id_col
            self$split_col_prefix <- split_col_prefix            
        },

        #' @description Read expression data into the matrix `expr_mat` and pheno data into the
        #' tibble `pheno_tbl`. Both will hold patients as rows.
        read = function()
        {  
            # extract values from self
            directory <- self$directory
            expr_file <- self$expr_file
            pheno_file <- self$pheno_file
            patient_id_col <- self$patient_id_col
            gene_id_col <- self$gene_id_col

            # read
            files <- c(expr_file, pheno_file)
            tbls <- list()
            for (i in 1:length(files)){
                full_path <- file.path(directory, files[i])
                tbls[[i]] <- readr::read_csv(full_path, show_col_types = FALSE)
            }
            expr_tbl <- tbls[[1]]
            pheno_tbl <- tbls[[2]]

            # check if identifier columns are there
            if(is.null(expr_tbl[[gene_id_col]])){
                stop("There is no column named ", gene_id_col, " in ", expr_file)
            }
            if(is.null(pheno_tbl[[patient_id_col]])){
                stop("There is no column named ", patient_id_col, " in ", pheno_file)
            }

            # expression to matrix
            gene_names <- tbls[["expr"]][[gene_id_col]]
            if(!elements_unique(gene_names)){
                stop("Column ", gene_id_col, " in ", expr_file, " holds duplicate entries.")
            }
            expr_mat <- expr_tbl |>
                dplyr::select(!dplyr::all_of(gene_id_col)) |>
                as.matrix() |> 
                t()
            colnames(expr_mat) <- gene_names

            # pheno: move patient ids into first column
            patient_ids <- pheno_tbl[[patient_id_col]]
            if(!elements_unique(patient_ids)){
                stop("Column ", patient_id_col, " in ", pheno_file, " holds duplicate entries.")
            }
            pheno_tbl <- pheno_tbl |>
                dplyr::relocate(dplyr::all_of(patient_id_col))

            self$expr_mat <- expr_mat
            self$pheno_tbl <- pheno_tbl
        }
    )
)

#' @title To every `ModelSpec` in a list of `ModelSpec`s, prepend a fixed directory to the
#' `directory` attribute
#' @description Often one specifies the models in general, for all data sets. If you fit the 
#' models to a specific data set (say `"mock"`), you might want to prepend a fixed directory 
#' like `"results/mock"` to the `directory` attribute of all `ModelSpec`s in the list.
#' @param model_list list of `ModelSpec`s.
#' @param dir string. The directory to prepend to the `directory` attribute of all 
#' `ModelSpec`s in `model_list`.
#' @return A list of `ModelSpec`s, with `dir` prepended to the `directory` attribute.
#' @export
prepend_to_directory <- function(
    model_list, 
    dir
){
    stopifnot(is.character(dir))
    for(i in seq_along(model_list)){
        model_list[[i]]$directory <- file.path(dir, model_list[[i]]$directory)
    }
    return(model_list)
}