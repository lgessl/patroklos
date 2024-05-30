data_prepare <- function(
    self,
    private,    
    model
){
    if(length(model$split_index) != 1)
        stop("Model must specify exactly one split index, but we have ", 
            length(model$split_index))
    if(length(model$time_cutoffs) != 1)
        stop("Model must specify exactly one time cutoff")
    if(is.null(self$cohort))
        stop("Cohort must be set to `'train'` or `'test'`")
    if(is.null(self$expr_mat) || is.null(self$pheno_tbl))
        stop("Data must have `expr_mat` and `pheno_tbl` read in")
    if(!all(rownames(self$expr_mat) == self$pheno_tbl[[self$patient_id_col]]))
        stop("Patient ids in expression and pheno data are not identical")

    x <- prepare_x(
        data = self,
        model = model
    )
    y <- prepare_y(
        data = self,
        model = model
    )
    y <- y[rownames(x), , drop = FALSE]

    return(list("x" = x, "y" = y))
}

# @title Generate the predictor matrix in a model-specific way
# @description Generate the numeric predictor matrix from the expression and
# (possibly) pheno data for a certain model
# @param expr_mat numeric matrix. The expression data, with patients as rows
# and genes as columns.
# @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
# columns.
# @param data Data S3 object. Specifications on the data (`expr_mat` and 
# `model`). See the the constructor `Data()` for details.
# @param model Model S3 object. Specifications on the model to prepare for. 
# See the the constructor `Model$new()` for details.
# @return A named numeric matrix with patients as rows and variables as columns. No 
# `NA`s in it.
prepare_x <- function(
    data,
    model
){
    # Extract
    patient_id_col <- data$patient_id_col
    include_from_continuous_pheno <- model$include_from_continuous_pheno
    include_from_discrete_pheno <- model$include_from_discrete_pheno
    x <- data$expr_mat
    pheno_tbl <- data$pheno_tbl
    # Check
    check_tbl_columns_exist(pheno_tbl, "pheno_tbl", c(include_from_discrete_pheno, 
        include_from_continuous_pheno))
    # Store levels of included discrete pheno variables for dichotomization
    pheno_tbl_ifdp <- pheno_tbl[, include_from_discrete_pheno, drop = FALSE]
    level_list <- lapply(pheno_tbl_ifdp, function(c) levels(as.factor(c)))

    # Subset pheno data ...
    # ... to cohort
    split_colname <- paste0(data$split_col_prefix, model$split_index)
    if(!split_colname %in% colnames(data$pheno_tbl))
        stop("Column ", split_colname, " not found in pheno table.")
    in_cohort_bool <- data$pheno_tbl[[split_colname]] == data$cohort
    if(all(in_cohort_bool) || all(!in_cohort_bool))
        stop("All patients are in the same cohort")
    x <- x[in_cohort_bool, , drop = FALSE]
    rnames_x <- rownames(x)
    pheno_tbl <- pheno_tbl[in_cohort_bool, , drop = FALSE]
    # ... to included variables
    pheno_tbl <- pheno_tbl[, c(include_from_continuous_pheno, 
        include_from_discrete_pheno), drop = FALSE]

    bind_continuous <- NULL
    bind_discrete <- NULL
    # Continuous pheno first
    if(!is.null(include_from_continuous_pheno)){
        bind_continuous <- as.matrix(pheno_tbl[, include_from_continuous_pheno, 
            drop = FALSE]) 
        if(!is.numeric(bind_continuous))
            stop("Variables specified in `include_from_continuous_pheno` must be numeric.")
        colnames(bind_continuous) <- paste0(colnames(bind_continuous), model$li_var_suffix) 
    }
    # Discrete pheno second
    if(!is.null(include_from_discrete_pheno)){
        pheno_tbl_ifdp <- pheno_tbl[, include_from_discrete_pheno, drop = FALSE]
        bind_discrete <- dichotomize_tibble(tbl = pheno_tbl_ifdp, level_list = 
            level_list)
        colnames(bind_discrete) <- paste0(colnames(bind_discrete), model$li_var_suffix)
    }

    # Combine into numeric matrix, the predictor matrix    
    if(!model$include_expr)
        x <- NULL
    x <- cbind(x, bind_continuous, bind_discrete)
    rownames(x) <- rnames_x

    # Impute
    if(!is.null(data$imputer)) {
        x_imp <- data$imputer(x)
        if (!all(colnames(x) == colnames(x_imp)) || 
            !all(rownames(x) == rownames(x_imp)))
            stop("Data$imputer() must return a matrix with the same dim names as ",
                "input matrix")
        if(!(all(x_imp[!is.na(x)] == x[!is.na(x)])))
            stop("Data$imputer() must not change non-NA values in input matrix")
        x <- x_imp
    }

    x <- x[stats::complete.cases(x), , drop = FALSE]
    attr(x, "li_var_suffix") <- model$li_var_suffix
    if(nrow(x) == 0){
        msg <- "No samples left in predictor matrix after "
        if (is.null(data$imputer))
            msg <- paste0(msg, "imputing and ")
        msg <- paste0(msg, "removing samples with NA.")
        stop(msg)
    }
    return(x)
}


# @title Generate the response matrix in a model-specific way
# @description Generate the numeric response matrix from the pheno data for a
# certain model
# @param pheno_tbl tibble. The pheno data, with patients as rows and variables as
# columns.
# @param data Data S3 object. Specifications on the data (`pheno_tbl`). See
# the constructor `Data()` for details.
# @param model Model S3 object. Specifications on the model to prepare for.
# See the constructor `Model$new()` for details.
# @return Named response matrix: a numeric matrix with patients as rows and variables 
# as columns.
# @details If `model$response_type == "binary"`, the response matrix will have one 
# column filled with 
# * `1` if progress is observed at a time <= `time_cutoff`,
# * `0` if progress or censoring is observed at a time > `time_cutoff`,
# * `NA` if censoring without progression is observed at a time <= `time_cutoff`.
prepare_y <- function(
    data,
    model
){

    # Extract
    time_to_event_col <- data$time_to_event_col
    event_col <- data$event_col
    patient_id_col <- data$patient_id_col
    time_cutoff <- model$time_cutoffs
    response_type <- model$response_type
    pheno_tbl <- data$pheno_tbl

    if(length(time_cutoff) > 1)
        stop("Can only handle one cutoff time.")

    if(response_type == "binary"){
        # flag patients censored before time_cutoff as NA
        na_bool <- (pheno_tbl[[time_to_event_col]] <= time_cutoff) & (pheno_tbl[[event_col]] == 0)
        y <- as.numeric(pheno_tbl[[time_to_event_col]] <= time_cutoff)
        dim(y) <- c(length(y), 1)
        rownames(y) <- pheno_tbl[[patient_id_col]]
        colnames(y) <- stringr::str_c("time_cutoff_", round(time_cutoff, 1))
        y[na_bool, ] <- NA
    } else if(response_type == "survival_censored"){
        y <- pheno_tbl[, c(time_to_event_col, event_col)] |> as.matrix()
        # Censor patients with time_to_event > time_cutoff at time_cutoff
        censor_bool <- y[, time_to_event_col] > time_cutoff
        y[censor_bool, 1] <- time_cutoff
        y[censor_bool, 2] <- 0
        rownames(y) <- pheno_tbl[[patient_id_col]]
        colnames(y) <- model$response_colnames
    }

    return(y)
}

mean_impute <- function(x) {
    x_imp <- sapply(seq(ncol(x)), function(j) {
        col <- x[, j]
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
    })
    colnames(x_imp) <- colnames(x)
    x_imp
}