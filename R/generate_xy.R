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
generate_predictor <- function(
    expr_mat,
    pheno_tbl,
    data,
    model
){
    # Extract
    patient_id_col <- data$patient_id_col
    include_from_continuous_pheno <- model$include_from_continuous_pheno
    include_from_discrete_pheno <- model$include_from_discrete_pheno

    x <- expr_mat
    patient_ids <- rownames(expr_mat) # Store for later
    bind_continuous <- NULL
    bind_discrete <- NULL

    check_consistent_patient_ids(
        stage = "before_generate_predictor", 
        expr = expr_mat, 
        pheno = pheno_tbl,
        data = data
    )

    # Continuous pheno first
    if(!is.null(include_from_continuous_pheno)){
        check_tbl_columns_exist(pheno_tbl, "pheno_tbl", include_from_continuous_pheno)
        bind_continuous <- pheno_tbl[, include_from_continuous_pheno, drop = FALSE] |> 
            as.matrix()
        if(!is.numeric(bind_continuous))
            stop("Variables specified in `include_from_continuous_pheno` must be numeric.")
        colnames(bind_continuous) <- colnames(bind_continuous) |> 
            stringr::str_c(model$li_var_suffix)
    }
    # Discrete pheno second
    if(!is.null(include_from_discrete_pheno)){
        check_tbl_columns_exist(pheno_tbl, "pheno_tbl", include_from_discrete_pheno)
        bind_discrete <- pheno_tbl[, include_from_discrete_pheno, drop = FALSE] |>
            tibble_to_binary()
        colnames(bind_discrete) <- colnames(bind_discrete) |>
            stringr::str_c(model$li_var_suffix)
    }

    # Combine into numeric matrix, the predictor matrix    
    x <- x |> cbind(bind_continuous, bind_discrete)
    rownames(x) <- patient_ids
    # Remove rows with NA (they are never useful)
    x <- x[stats::complete.cases(x), , drop = FALSE]
    attr(x, "li_var_suffix") <- model$li_var_suffix

    if(nrow(x) == 0)
        stop("No samples left in predictor matrix after subsetting to cohort and removing NAs")

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
generate_response <- function(
    pheno_tbl,
    data,
    model
){

    # Extract
    time_to_event_col <- data$time_to_event_col
    event_col <- data$event_col
    patient_id_col <- data$patient_id_col
    time_cutoff <- model$time_cutoffs
    response_type <- model$response_type

    if(length(time_cutoff) > 1)
        stop("Can only handle one cutoff time.")

    if(response_type == "binary"){
        # flag patients censored before time_cutoff as NA
        na_bool <- (pheno_tbl[[time_to_event_col]] <= time_cutoff) & (pheno_tbl[[event_col]] == 0)
        y <- pheno_tbl[[time_to_event_col]] <= time_cutoff
        y <- as.numeric(y)
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