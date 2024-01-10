# (Internally) construct a ModelSpec S3 object
new_ModelSpec <- function(
    name,
    fitter,
    cutoff_times,
    split_index,
    optional_fitter_args,
    response_type,
    response_colnames,
    include_from_continuous_pheno,
    include_from_discrete_pheno,
    append_to_includes,
    pheno_regexp,
    base_dir,
    directory,
    create_directory,
    plot_fname,
    plot_ncols,
    fit_fname
){
    stopifnot(is.character(name))
    check_fitter(fitter, optional_fitter_args)
    stopifnot(is.numeric(cutoff_times))
    stopifnot(is.numeric(split_index))
    stopifnot(is.character(response_type))
    stopifnot(is.character(response_colnames))
    stopifnot(is.character(include_from_continuous_pheno) || is.null(include_from_continuous_pheno))
    stopifnot(is.character(include_from_discrete_pheno) || is.null(include_from_discrete_pheno))
    stopifnot(is.character(append_to_includes))
    stopifnot(is.character(base_dir))
    stopifnot(is.character(directory))
    stopifnot(is.logical(create_directory))
    stopifnot(is.numeric(cutoff_times) || is.null(cutoff_times))
    stopifnot(is.character(plot_fname))
    stopifnot(is.numeric(plot_ncols))
    stopifnot(is.character(fit_fname))

    model_spec_list <- list(
        "name" = name,
        "fitter" = fitter,
        "cutoff_times" = cutoff_times,
        "split_index" = split_index,
        "optional_fitter_args" = optional_fitter_args,
        "response_type" = response_type,
        "response_colnames" = response_colnames,
        "include_from_continuous_pheno" = include_from_continuous_pheno,
        "include_from_discrete_pheno" = include_from_discrete_pheno,
        "append_to_includes" = append_to_includes,
        "cutoff_times" = cutoff_times,
        "base_dir" = base_dir,
        "directory" = directory,
        "create_directory" = create_directory,
        "plot_fname" = plot_fname,
        "plot_ncols" = plot_ncols,
        "fit_fname" = fit_fname
    )
    return(structure(model_spec_list, class = "ModelSpec"))
}

#' @title Construct a ModelSpec S3 object
#' @description A ModelSpec object holds all the tools and information needed to fit and
#' store models. Most importantly, it is passed as an argument to [`training_camp()`]. 
#' Its base object is a list.
#' @param name string. A telling name for the model.
#' @param fitter function. The model fitting function to be used. Must take `x` and
#' `y` as first two positional arguments. Further arguments can be passed via
#' `optional_fitter_args` (below). Its return value must be an S3 object with a `plot()` 
#' method, and (ideally, for assessment) with a `predict()` method. Default is `NULL`.
#' @param cutoff_times numeric vector.
#' * If `response_type == "survival_censored"`: For every value in `cutoff_times`, censor
#' all patients where the event ouccured after `cutoff_times` at this value and train the 
#' specified model.
#' * If `response_type == "binary"`: For every value in `cutoff_times`, binarize the 
#' outcome depending on whether it occured before or after this value, and train the 
#' specified model.
#' @param split_index integer vector. Split the given data into training and test samples 
#' `length(split_index)` times, i.e., every index in `split_index` will get its own split. 
#' @param optional_fitter_args list. Optional arguments passed to `fitter`, e.g. alpha 
#' in case of an elastic net. Default is `list()`, i.e., no arguments other than `x`, `y`
#' passed to `fitter`.
#' @param response_type string. The type of response to be used. One of `"binary"` or
#' `"survival_censored"`. Default is `NULL`.
#' @param response_colnames string vector of length 2. If `response_type == "survival_censored"`,
#' use as column names for the response matrix. 
#' * The first element is the name of the column holding the time until the event or 
#' censoring, and 
#' * the second one is the anme of the column holding the event status (1 = event, 0 =
#' censoring). 
#' Default is `c("time", "status")`.
#' @param include_from_continuous_pheno vector of strings. The names of the 
#' *continuous* variables in the pheno data (to be) included in the predictor matrix. The
#' values will be coerced to numeric. Default is `NULL`, which means no continuous pheno
#' variables are or will be included.
#' @param include_from_discrete_pheno vector of strings. The names of the *discrete*
#' variables in the pheno data (to be) included in the predictor matrix. A discrete
#' variable with n levels will be converted to n-1 binary variables. Default is `NULL`,
#' which means no discrete pheno variables are or will be included.
#' @param append_to_includes string. Append this to the names of features from the pheno
#' data when adding them to the predictor matrix. Default is `"++"`.
#' @param base_dir string. The base directory to store the model in. See `directory` below
#' on how it is used to automatically set `directory`. Default is `"."`.
#' @param directory string. The directory to store the models in. For every value in 
#' `cutoff_times`, find the corresponding model in a subdirectory named after this value. 
#' Default is `NULL`, in which case is is set to `file.path(base_dir, name)`.
#' @param create_directory logical. Whether to create `directory` if it does not exist, yet. 
#' Default is `TRUE`.
#' @param plot_fname string. Store the plot resulting from `plot(fit_obj)` in `directory`
#' under this name. Default is `"training_error.pdf"`.
#' @param plot_ncols integer. The number of columns in the plot. Default is `2`.
#' @param fit_fname string. The name of the model-fits file inside `directory`.
#' Default is `"fit_obj.rds"`.
#' @return A ModelSpec S3 object.
#' @details Strictly speaking, one `ModelSpec` instance holds the instructions to fit
#' `length(cutoff_times) * length(split_index)` models. In terms of storing and assessing models,
#' we consider the models obtained via repeated splitting according to `split_index` as one 
#' model; repeated splitting serves the purpose of getting more reliable estimates of its
#' performance. We view models obtained via different values of `cutoff_times`, in 
#' contrast, as different models; e.g., we can compare them against one another in an 
#' assessment.
#' @export
ModelSpec <- function(
    name,
    directory,
    fitter,
    split_index,
    cutoff_times,
    optional_fitter_args = NULL,
    response_type = c("binary", "survival_censored"),
    response_colnames = c("time", "status"),
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    append_to_includes = "++",
    base_dir = ".",
    create_directory = TRUE,
    plot_fname = "training_error.pdf",
    plot_ncols = 2,
    fit_fname = "fit_obj.rds"
){
    if(is.null(directory)){
        directory <- file.path(base_dir, name)
    }
    response_type <- match.arg(response_type)
    stopifnot(all(cutoff_times >= 0))
    stopifnot(split_index >= 1)
    model_spec <- new_ModelSpec(
        name = name,
        fitter = fitter,
        cutoff_times = cutoff_times,
        split_index = split_index,
        optional_fitter_args = optional_fitter_args,
        response_type = response_type,
        response_colnames = response_colnames,
        include_from_continuous_pheno = include_from_continuous_pheno,
        include_from_discrete_pheno = include_from_discrete_pheno,
        append_to_includes = append_to_includes,
        directory = directory,
        base_dir = base_dir,
        create_directory = create_directory,
        plot_fname = plot_fname,
        plot_ncols = plot_ncols,
        fit_fname = fit_fname
    )
    return(model_spec)
}