model_initialize <- function(self, private, name, fitter, directory, split_index, 
    time_cutoffs, hyperparams, response_type, response_colnames, 
    include_from_continuous_pheno, include_from_discrete_pheno, include_expr, 
    li_var_suffix, create_directory, plot_file, plot_ncols,
    plot_title_line, fit_file, continuous_output){

    response_type <- match.arg(response_type, c("binary", "survival_censored"))
    stopifnot(all(time_cutoffs >= 0))
    stopifnot(split_index >= 1)
    stopifnot(is.character(name))
    check_fitter(fitter, hyperparams)
    stopifnot(is.numeric(time_cutoffs))
    stopifnot(is.numeric(split_index))
    stopifnot(is.character(response_type))
    stopifnot(is.character(response_colnames))
    stopifnot(is.character(include_from_continuous_pheno) || is.null(include_from_continuous_pheno))
    stopifnot(is.character(include_from_discrete_pheno) || is.null(include_from_discrete_pheno))
    stopifnot(is.logical(include_expr))
    stopifnot(is.character(li_var_suffix))
    stopifnot(is.character(directory))
    stopifnot(is.logical(create_directory))
    stopifnot(is.numeric(time_cutoffs) || is.null(time_cutoffs))
    stopifnot(is.character(plot_file))
    stopifnot(is.numeric(plot_ncols))
    stopifnot(is.numeric(plot_title_line) || is.null(plot_title_line))
    stopifnot(is.character(fit_file))
    stopifnot(is.logical(continuous_output) || is.null(continuous_output))

    self$name <- name
    self$directory <- directory
    self$fitter <- fitter
    self$split_index <- split_index
    self$time_cutoffs <- time_cutoffs
    self$hyperparams <- hyperparams
    self$response_type <- response_type
    self$response_colnames <- response_colnames
    self$include_from_continuous_pheno <- include_from_continuous_pheno
    self$include_from_discrete_pheno <- include_from_discrete_pheno
    self$include_expr <- include_expr
    self$li_var_suffix <- li_var_suffix
    self$create_directory <- create_directory
    self$plot_file <- plot_file
    self$plot_ncols <- plot_ncols
    self$plot_title_line <- plot_title_line
    self$fit_file <- fit_file
    self$fits <- vector("list", length(split_index))
    self$continuous_output <- continuous_output

    invisible(self)
}

model_at_time_cutoff <- function(
    self,
    private,
    time_cutoff
){
    stopifnot(is.numeric(time_cutoff))
    if(length(time_cutoff) != 1) stop("time_cutoff must be a single number")
    if(!time_cutoff %in% self$time_cutoffs){
        stop("time_cutoff must be one of model$time_cutoffs")
    }
    model_cutoff <- self$clone()
    model_cutoff$name <- paste0(self$name, "@", time_cutoff)
    model_cutoff$time_cutoffs <- time_cutoff
    model_cutoff$directory <- file.path(
        self$directory, 
        stringr::str_replace(as.character(time_cutoff), "\\.", "-")
    )
    invisible(model_cutoff)
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
    invisible(NULL)
}