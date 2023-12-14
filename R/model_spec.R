# (Internally) construct a ModelSpec S3 object
new_ModelSpec <- function(
    name,
    fitter,
    optional_fitter_args,
    response_type,
    include_from_continuous_pheno,
    include_from_discrete_pheno,
    append_to_includes,
    pheno_regexp,
    base_dir,
    save_dir,
    create_save_dir,
    pfs_leq,
    plot_fname,
    fit_fname
){
    check_fitter(fitter, optional_fitter_args)
    response_type <- match.arg(response_type, c("binary", "survival_censored"))
    stopifnot(is.character(name))
    stopifnot(is.character(include_from_continuous_pheno) || is.null(include_from_continuous_pheno))
    stopifnot(is.character(include_from_discrete_pheno) || is.null(include_from_discrete_pheno))
    stopifnot(is.character(append_to_includes))
    stopifnot(is.character(base_dir))
    stopifnot(is.character(save_dir))
    stopifnot(is.logical(create_save_dir))
    stopifnot(is.numeric(pfs_leq) || is.null(pfs_leq))
    stopifnot(is.character(plot_fname))
    stopifnot(is.character(fit_fname))

    model_spec_list <- list(
        "name" = name,
        "fitter" = fitter,
        "optional_fitter_args" = optional_fitter_args,
        "response_type" = response_type,
        "include_from_continuous_pheno" = include_from_continuous_pheno,
        "include_from_discrete_pheno" = include_from_discrete_pheno,
        "append_to_includes" = append_to_includes,
        "pfs_leq" = pfs_leq,
        "base_dir" = base_dir,
        "save_dir" = save_dir,
        "create_save_dir" = create_save_dir,
        "plot_fname" = plot_fname,
        "fit_fname" = fit_fname
    )
    return(structure(model_spec_list, class = "ModelSpec"))
}

#' @title Construct a ModelSpec S3 object
#' @description A ModelSpec object holds all the tools and information needed to fit and
#' store a model. It is used as an argument to `fit()` and `prepare_and_fit()`. Its base
#' object is a list containing:
#' @param name string. A telling name for the model.
#' @param fitter function. The model fitting function to be used. Must take `x` and
#' `y` as first two positional arguments. Further arguments can be passed via
#' `optional_fitter_args` (below). Its return value must be an S3 object with a `plot()` 
#' method, and (ideally, for assessment) with a `predict()` method. Default is `NULL`.
#' @param optional_fitter_args list. Optional arguments passed to `fitter`, e.g. alpha 
#' in case of an elastic net. Default is `list()`, i.e., no arguments other than `x`, `y`
#' passed to `fitter`.
#' @param response_type string. The type of response to be used. One of `"binary"` or
#' `"survival_censored"`. Default is `NULL`.
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
#' @param pfs_leq numeric. Only used if `response_type == "binary"`. The value of
#' progression-free survival (PFS) below which samples are considered high-risk. Default
#' is `2.0`.
#' @param base_dir string. The base directory to store the model in. See `save_dir` below
#' on how it is used to automatically set `save_dir`. Default is `"."`.
#' @param save_dir string. The directory in which to store the model in. Default is `NULL`, 
#' in which case is is set to `file.path(base_dir, name)`.
#' @param create_save_dir logical. Whether to create `save_dir` if it does not exist, yet. 
#' Default is `TRUE`.
#' @param plot_fname string. Store the plot resulting from `plot(fit_obj)` in `save_dir`
#' under this name. Default is `"training_error.pdf"`.
#' @param fit_fname string. The name of the model-fit file inside `save_dir`.
#' Default is `"fit_obj.rds"`.
#' @return A ModelSpec S3 object.
#' @export
ModelSpec <- function(
    name,
    fitter,
    optional_fitter_args = NULL,
    response_type = NULL,
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    append_to_includes = "++",
    pfs_leq = 2.0,
    base_dir = ".",
    save_dir = NULL,
    create_save_dir = TRUE,
    plot_fname = "training_error.pdf",
    fit_fname = "fit_obj.rds"
){
    if(is.null(save_dir)){
        save_dir <- file.path(base_dir, name)
    }
    model_spec <- new_ModelSpec(
        name = name,
        fitter = fitter,
        optional_fitter_args = optional_fitter_args,
        response_type = response_type,
        include_from_continuous_pheno = include_from_continuous_pheno,
        include_from_discrete_pheno = include_from_discrete_pheno,
        append_to_includes = append_to_includes,
        save_dir = save_dir,
        base_dir = base_dir,
        create_save_dir = create_save_dir,
        pfs_leq = pfs_leq,
        plot_fname = plot_fname,
        fit_fname = fit_fname
    )
    return(model_spec)
}