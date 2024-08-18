model_initialize <- function(self, private, name, fitter, directory, 
    time_cutoffs, val_error_fun, hyperparams, 
    include_from_continuous_pheno, include_from_discrete_pheno, include_expr, 
    li_var_suffix, create_directory, file, 
    combine_n_max_categorical_features, combined_feature_min_positive_ratio, enable_imputation){

    stopifnot(all(time_cutoffs >= 0))
    stopifnot(is.character(name))
    check_fitter(fitter, hyperparams)
    stopifnot(is.numeric(time_cutoffs))
    stopifnot(is.function(val_error_fun))
    stopifnot(is.character(include_from_continuous_pheno) || is.null(include_from_continuous_pheno))
    stopifnot(is.character(include_from_discrete_pheno) || is.null(include_from_discrete_pheno))
    stopifnot(is.logical(include_expr))
    stopifnot(is.character(li_var_suffix))
    stopifnot(is.character(directory))
    stopifnot(is.logical(create_directory))
    stopifnot(is.numeric(time_cutoffs) || is.null(time_cutoffs))
    stopifnot(is.character(file))
    combine_n_max_categorical_features <- as.integer(round(combine_n_max_categorical_features))
    if (any(combine_n_max_categorical_features < 1))
        stop("combine_n_max_categorical_features must be at least 1.")
    stopifnot(is.numeric(combined_feature_min_positive_ratio) && 
        combined_feature_min_positive_ratio >= 0 && 
        combined_feature_min_positive_ratio <= 1)
    stopifnot(is.logical(enable_imputation))

    self$name <- name
    self$directory <- directory
    self$fitter <- fitter
    self$time_cutoffs <- time_cutoffs
    self$val_error_fun <- val_error_fun
    self$hyperparams <- hyperparams
    self$include_from_continuous_pheno <- include_from_continuous_pheno
    self$include_from_discrete_pheno <- include_from_discrete_pheno
    self$include_expr <- include_expr
    self$li_var_suffix <- li_var_suffix
    self$create_directory <- create_directory
    self$file <- file
    self$combine_n_max_categorical_features <- combine_n_max_categorical_features
    self$combined_feature_min_positive_ratio <- combined_feature_min_positive_ratio
    self$enable_imputation <- enable_imputation

    invisible(self)
}

#' @title To every `Model` in a list of `Model`s, prepend a fixed directory to the
#' `directory` attribute
#' @description This function supports the paradigm of specifying the `Model`s for 
#' for *all* data sets at one place and them quickly adopting for a specific data 
#' set. We modify the Model objects *in place*.
#' @param model_list list of `Model`s.
#' @param prefix character Prepend this directory to the `directory` attribute of each
#' `Model` in `model_list`.
#' @export
prepend_to_directory <- function(
    model_list, 
    prefix
){
    stopifnot(is.character(prefix))
    lapply(model_list, function(model) {
        stopifnot(inherits(model, "Model"))
        model$directory <- file.path(prefix, model$directory)
    })
    n_unique_dir <- length(unique(sapply(model_list, function(x) x$directory)))
    if (n_unique_dir != length(model_list)) 
        stop("Model directories are not unique.")
    n_unique_names <- length(unique(sapply(model_list, function(x) x$name)))
    if (n_unique_names != length(model_list)) 
        stop("Model names are not unique.")
    invisible(NULL)
}