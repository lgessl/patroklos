model_initialize <- function(self, private, name, fitter, directory, split_index, 
    time_cutoffs, optional_fitter_args, response_type, response_colnames, 
    include_from_continuous_pheno, include_from_discrete_pheno, 
    append_to_includes, create_directory, plot_file, plot_ncols,
    plot_title_line, fit_file){

    response_type <- match.arg(response_type, c("binary", "survival_censored"))
    stopifnot(all(time_cutoffs >= 0))
    stopifnot(split_index >= 1)
    stopifnot(is.character(name))
    check_fitter(fitter, optional_fitter_args)
    stopifnot(is.numeric(time_cutoffs))
    stopifnot(is.numeric(split_index))
    stopifnot(is.character(response_type))
    stopifnot(is.character(response_colnames))
    stopifnot(is.character(include_from_continuous_pheno) || is.null(include_from_continuous_pheno))
    stopifnot(is.character(include_from_discrete_pheno) || is.null(include_from_discrete_pheno))
    stopifnot(is.character(append_to_includes))
    stopifnot(is.character(directory))
    stopifnot(is.logical(create_directory))
    stopifnot(is.numeric(time_cutoffs) || is.null(time_cutoffs))
    stopifnot(is.character(plot_file))
    stopifnot(is.numeric(plot_ncols))
    stopifnot(is.numeric(plot_title_line) || is.null(plot_title_line))
    stopifnot(is.character(fit_file))

    self$name <- name
    self$directory <- directory
    self$fitter <- fitter
    self$split_index <- split_index
    self$time_cutoffs <- time_cutoffs
    self$optional_fitter_args <- optional_fitter_args
    self$response_type <- response_type
    self$response_colnames <- response_colnames
    self$include_from_continuous_pheno <- include_from_continuous_pheno
    self$include_from_discrete_pheno <- include_from_discrete_pheno
    self$append_to_includes <- append_to_includes
    self$create_directory <- create_directory
    self$plot_file <- plot_file
    self$plot_ncols <- plot_ncols
    self$plot_title_line <- plot_title_line
    self$fit_file <- fit_file
    self$fits <- vector("list", length(split_index))

    invisible(self)
}

model_at_time_cutoff <- function(
    self,
    private,
    time_cutoff
){
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
#' @description Often one specifies the models in general, for all data sets. If you fit the 
#' models to a specific data set (say `"mock"`), you might want to prepend a fixed directory 
#' like `"results/mock"` to the `directory` attribute of all `Model`s in the list.
#' @param model_list list of `Model`s.
#' @param dir string. The directory to prepend to the `directory` attribute of all 
#' `Model`s in `model_list`.
#' @return A list of `Model`s, with `dir` prepended to the `directory` attribute.
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