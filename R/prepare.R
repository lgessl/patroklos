data_prepare <- function(
    self,
    private,    
    model,
    quiet
){
    if(length(model$split_index) != 1)
        stop("Model must specify exactly one split index, but we have ", 
            length(model$split_index))
    if(!is.character(self$cohort))
        stop("data$cohort must be set")
    if(is.null(self$expr_mat) || is.null(self$pheno_tbl))
        stop("Data must have `expr_mat` and `pheno_tbl` read in")
    if(!all(rownames(self$expr_mat) == self$pheno_tbl[[self$patient_id_col]]))
        stop("Patient ids in expression and pheno data are not identical")

    x <- prepare_x(
        data = self,
        model = model,
        quiet = quiet
    )
    y_cox <- as.matrix(self$pheno_tbl[, c(self$time_to_event_col, self$event_col)])
    rownames(y_cox) <- self$pheno_tbl[[self$patient_id_col]]
    colnames(y_cox) <- model$response_colnames

    # Subset y_cox a bit: we don't need outcomes we can't predict for
    y_cox <- y_cox[rownames(x), , drop = FALSE]

    list(x = x, y_cox = y_cox)
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
    model,
    quiet
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

    bind_continuous <- x[, 0, drop = FALSE]
    bind_discrete <- x[, 0, drop = FALSE]
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
        bind_discrete <- dichotomize_tibble(tbl = pheno_tbl_ifdp)
    }

    # Combine into numeric matrix (the first time: for imputation)
    rnames_x <- rownames(x)
    if(!model$include_expr)
        x <- x[, 0, drop = FALSE]
    discrete_col_bool <- c(logical(ncol(x) + ncol(bind_continuous)), 
        !logical(ncol(bind_discrete)))
    x <- cbind(x, bind_continuous, bind_discrete)
    rownames(x) <- rnames_x

    # Prepare subsetting to cohort
    split_colname <- paste0(data$split_col_prefix, model$split_index)
    if(!split_colname %in% colnames(data$pheno_tbl))
        stop("Column ", split_colname, " not found in pheno table.")
    in_cohort_bool <- stringr::str_detect(data$pheno_tbl[[split_colname]], data$cohort)
    if(all(in_cohort_bool) || all(!in_cohort_bool))
        stop("All patients are in the same cohort")

    # Impute
    if(!is.null(data$imputer)) {
        x_imp <- x
        cohort_col <- data$pheno_tbl[[split_colname]]
        for (cohort in unique(cohort_col)) {
            take_bool <- cohort_col == cohort
            x_imp[take_bool, ] <- data$imputer(x[take_bool, , drop = FALSE])
        }
        if (!all(colnames(x) == colnames(x_imp)) || 
            !all(rownames(x) == rownames(x_imp)))
            stop("Data$imputer() must return a matrix with the same dim names as ",
                "input matrix")
        if(!(all(x_imp[!is.na(x)] == x[!is.na(x)])))
            stop("Data$imputer() must not change non-NA values in input matrix")
        x <- x_imp
    }

    # Combine discrete features
    if (!is.null(include_from_discrete_pheno)) {
        bind_discrete_wide <- combine_features(
                x = x[, discrete_col_bool, drop = FALSE], 
                combine_n_max_features = max(model$combine_n_max_categorical_features),
                combined_feature_positive_ratio = model$combined_feature_min_positive_ratio, 
                original_cnames = include_from_discrete_pheno
        )
        if (!quiet)
            message("** Including ", ncol(bind_discrete_wide), 
                " combined categorical features.")
        colnames(bind_discrete_wide) <- paste0(colnames(bind_discrete_wide),
            model$li_var_suffix)
        x <- cbind(x[, !discrete_col_bool, drop = FALSE], bind_discrete_wide)
        rownames(x) <- rnames_x
    } 
    
    # Actually and finally subset to cohort
    x <- x[in_cohort_bool, , drop = FALSE]
    if (!quiet)
        message("** Dimensions of prepared x: (", paste(dim(x), collapse = ", "), ").")
    attr(x, "li_var_suffix") <- model$li_var_suffix
    return(x)
}

# Turn y_cox into binary response y_bin according to time_cutoff
binarize_y <- function(y_cox, time_cutoff, pivot_time_cutoff) {
    # Flag patients censored before time_cutoff as NA
    if (is.infinite(time_cutoff)) time_cutoff <- pivot_time_cutoff
    na_bool <- (y_cox[, 1] <= time_cutoff) & (y_cox[, 2] == 0)
    y_bin <- as.matrix((y_cox[, 1] <= time_cutoff) * 1)
    y_bin[na_bool, ] <- NA
    y_bin
}

# In y_cox, censor samples with time to event > time_cutoff at time cutoff
censor_y <- function(y_cox, time_cutoff) {
    censor_bool <- y_cox[, 1] > time_cutoff
    y_cox[censor_bool, 1] <- time_cutoff
    y_cox[censor_bool, 2] <- 0
    y_cox
}

#' @title Impute missing values in a matrix by column means
#' @description For categorical variables (which we expect to be dichotomized 
#' already), we also simply compute the mean, which is the (marginal) Bernoullie 
#' probability to be in the positive class (1).
#' @param x numeric matrix. Sample as rows, variables as columns. Impute for 
#' `NA`s in the columns.
#' @export
mean_impute <- function(x) {
    x_imp <- sapply(seq(ncol(x)), function(j) {
        col <- x[, j]
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
    })
    colnames(x_imp) <- colnames(x)
    x_imp
}

combine_features <- function(x, combine_n_max_features, 
    combined_feature_positive_ratio, original_cnames) {

    # Get list of all possible combinations of columns
    cand_comb <- replicate(combine_n_max_features, seq(ncol(x)), simplify = FALSE)
    cand_comb <- expand.grid(cand_comb)
    cand_comb <- apply(cand_comb, 1, function(r) sort(unique(r)), simplify = FALSE)
    cand_comb <- unique(cand_comb)

    # No intra-feature combinations
    comb_cnames_list <- lapply(cand_comb, function(c) colnames(x)[c])
    only_inter_feat <- sapply(comb_cnames_list, function(cc) 
        all(sapply(original_cnames, function(oc) sum(stringr::str_detect(cc, oc)) <= 1)))
    cand_comb <- cand_comb[only_inter_feat]
    comb_cnames_list <- comb_cnames_list[only_inter_feat]
    comb_cnames <- sapply(comb_cnames_list, function(c) paste(c, collapse = "&"))

    # Muliplication equates to AND
    x_wide <- sapply(cand_comb, function(comb) apply(x[, comb, drop = FALSE], 1, 
        prod))
    x_wide <- as.matrix(x_wide)
    colnames(x_wide) <- comb_cnames

    # Only keep columns with sufficient positive ratio
    x_wide[, colMeans(x_wide, na.rm = TRUE) >= combined_feature_positive_ratio, 
        drop = FALSE]
}

# Remove columns with combinations of comprising more than
# combine_n_max_categorical_features
trim_combos <- function(x, combine_n_max_categorical_features) {
    keep <- sapply(stringr::str_split(colnames(x), "&"), length) <= 
        combine_n_max_categorical_features
    x_slim <- x[, keep]
    attr(x_slim, "li_var_suffix") <- attr(x, "li_var_suffix")
    x_slim
}