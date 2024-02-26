#' @title Quality control before handing over the predictor and response 
#' matrix over to the model fitting function
#' @description Check types of predictor and response matrix, their 
#' consistency and availability.
#' @param x A numeric matrix with dimnames holding predictors (as rows).
#' @param y A numeric matrix with dimnames holding responses (as rows).
#' @details Intended to be used after `generate_predictor()` and 
#' `generate_response()`. `prepare()` and, in particular, `prepare_and_fit()`
#' call this function.
qc_prefit <- function(
    x,
    y
){
    # check class and type of x and y
    x_y_list <- list("x" = x, "y" = y)
    for(x_y in names(x_y_list)){
        obj <- x_y_list[[x_y]]
        if(!is.matrix(obj))
            stop(
                x_y, " must be a matrix. But its S3 class is (", 
                stringr::str_c(sloop::s3_class(obj), collapse = ", "), ")"
            )
        if(!is.numeric(obj))
            stop(
                x_y, " must be a numeric. But its S3 class is (", 
                stringr::str_c(sloop::s3_class(obj), collapse = ", "), ")"
            )
    }
    check_consistent_patient_ids(
        stage = "after_generate_xy",
        expr = x,
        pheno = y
    )
    check_available(
        x = x,
        y = y
    )
}


#' @title Quality control at the end of preprocessing step
#' @description Check if the expression and pheno tibble have the format
#' specified in the `Data` and match to one another.
#' @param expr_tbl A tibble holding the expression data (see `Data()`
#' for details).
#' @param pheno_tbl A tibble holding the pheno data (see `Data()`
#' for details).
#' @param data A `Data` object referring to `expr_tbl` and `pheno_tbl`.
#' @param check_default logical. If `TRUE`, check if all default arguments of
#' `Data()` apply (i.e., replace `data` with a `Data()` except for 
#' the `directory` attribute). Default is `FALSE`, in which case we use `data`.
#' @details Intended to be used at the end of preprocessing and before applying
#' `split_data()`.
#' @export
qc_preprocess <- function(
    expr_tbl,
    pheno_tbl,
    data,
    check_default = FALSE
){
    if(check_default){
        directory <- data$directory
        data <- Data(
            name = "default",
            directory = directory,
            train_prop = .7
        )
    }

    # Extract
    directory <- data$directory
    expr_file <- data$expr_file
    pheno_file <- data$pheno_file
    gene_id_col <- data$gene_id_col
    patient_id_col <- data$patient_id_col
    time_to_event_col <- data$time_to_event_col
    event_col <- data$event_col
    benchmark_col <- data$benchmark_col

    # Check if files exist
    expr_file <- file.path(directory, expr_file)
    pheno_file <- file.path(directory, pheno_file)
    if(!file.exists(expr_file)){
        stop("Expression file ", expr_file, " does not exist.")
    }
    if(!file.exists(pheno_file)){
        stop("Pheno file ", pheno_file, " does not exist.")
    }

    # Expression tibble
    if(names(expr_tbl)[1] != gene_id_col){
        stop("First column of expression tibble must be ", gene_id_col)
    }
    if(!is.character(expr_tbl[[gene_id_col]])){
        stop("Gene ids must be characters.")
    }
    if(!all(sapply(expr_tbl[, -1], is.numeric))){
        stop("Expression values must be numeric.")
    }
    if(any(is.na(expr_tbl[, -1]))){
        stop("Expression values must not contain missing values.")
    }

    # Pheno tibble
    check_tbl_columns_exist(
        tbl = pheno_tbl,
        tbl_name = "pheno_tbl",
        col_names = c(
            patient_id_col,
            time_to_event_col,
            event_col,
            benchmark_col
        )
    )
    if(names(pheno_tbl)[1] != patient_id_col){
        stop("First column of pheno tibble must be ", patient_id_col)
    }
    if(!is.character(pheno_tbl[[patient_id_col]]) || !elements_unique(pheno_tbl[[patient_id_col]])){
        stop("Patient ids must be unique characters.")
    }
    if(!is.numeric(pheno_tbl[[time_to_event_col]])){
        stop("PFS values must be numeric and not contain missing values.")
    }
    if(any(is.na(pheno_tbl[[time_to_event_col]]))){
        warning("PFS values contain missing values.")
    }
    if(!is.numeric(pheno_tbl[[event_col]]) || 
        !all(pheno_tbl[[event_col]] %in% c(0, 1))){
        stop("Progression values must be numeric and either 1 (progression) or 0 (no progression).")
    }
    if(!is.null(benchmark_col) && 
        stringr::str_detect(benchmark_col, stringr::regex("ipi", ignore_case = TRUE))){
        if(!is.numeric(pheno_tbl[[benchmark_col]]) || !all(pheno_tbl[[benchmark_col]] %in% c(0:5, NA))){
            stop("IPI values must be numeric values in 0:5")
        }
    }

    # Both
    check_consistent_patient_ids(
        stage = "preprocessing",
        expr = expr_tbl,
        pheno = pheno_tbl,
        data = data
    )
}