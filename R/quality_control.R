# @title Quality control before handing over the predictor and response 
# matrix over to the model fitting function
# @description Check types of predictor and response matrix, their 
# consistency and availability.
# @param x A numeric matrix with dimnames holding predictors (as rows).
# @param y A numeric matrix with dimnames holding responses (as rows).
# @details Intended to be used after `prepare_x()` and 
# `prepare_y()`.
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
    if (is.null(colnames(x))) {
        stop("x has no column names")
    }
    if (is.null(attr(x, "li_var_suffix"))) {
        stop("x has no li_var_suffix attribute")
    }
}


data_qc_preprocess <- function(self, private, expr_tbl){

    if(is.null(self$pheno_tbl))
        stop("pheno_tbl attribute of Data object must not be NULL")

    # Expression tibble
    if(names(expr_tbl)[1] != self$gene_id_col)
        stop("First column of expression tibble must be ", self$gene_id_col)
    if(!is.character(expr_tbl[[self$gene_id_col]]) || 
        !elements_unique(expr_tbl[[self$gene_id_col]]))
        stop("Gene ids must be unique characters.")
    if(!all(sapply(expr_tbl[, -1], is.numeric)))
        stop("Expression values must be numeric.")
    if(any(is.na(expr_tbl[, -1])))
        stop("Expression values must not contain missing values")

    # Pheno tibble
    check_tbl_columns_exist(
        tbl = self$pheno_tbl,
        tbl_name = "pheno_tbl",
        col_names = c(
            self$patient_id_col,
            self$time_to_event_col,
            self$event_col,
            self$benchmark_col
        )
    )
    if(names(self$pheno_tbl)[1] != self$patient_id_col)
        stop("First column of pheno tibble must be ", self$patient_id_col)
    if(is.numeric(self$pheno_tbl[[self$patient_id_col]]))
        self$pheno_tbl[[self$patient_id_col]] <- 
            as.character(self$pheno_tbl[[self$patient_id_col]])
    if(!is.character(self$pheno_tbl[[self$patient_id_col]]) || 
        !elements_unique(self$pheno_tbl[[self$patient_id_col]]))
        stop("Patient ids must be unique characters")
    if(!is.numeric(self$pheno_tbl[[self$time_to_event_col]]))
        stop("Time-to-event values must be numeric")
    if(any(is.na(self$pheno_tbl[[self$time_to_event_col]])))
        warning("Time-to-event values contain missing values")
    if(!is.numeric(self$pheno_tbl[[self$event_col]]) || 
        !all(self$pheno_tbl[[self$event_col]] %in% c(0, 1)))
        stop("Progression values must be numeric and either 1 (progression) or 0 (no progression).")
    if(!is.null(self$benchmark_col) && 
        stringr::str_detect(self$benchmark_col, stringr::regex("ipi", ignore_case = TRUE))){
        if(!is.numeric(self$pheno_tbl[[self$benchmark_col]]) || 
            !all(self$pheno_tbl[[self$benchmark_col]] %in% c(0:5, NA))){
            stop("IPI values must be numeric values in 0:5")
        }
    }

    # Both
    check_consistent_patient_ids(
        stage = "preprocessing",
        expr = expr_tbl,
        pheno = self$pheno_tbl,
        data = self
    )
}