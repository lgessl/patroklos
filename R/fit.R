#' @title Fit and store a model
#' @description Given a model and data in a model-specific form, fit this model 
#' and store it
#' @param x numeric matrix. The predictor matrix with samples in rows and their names
#' as row names.
#' @param y numeric matrix. The response matrix with samples in rows and their names
#' as row names.
#' @param fitter function. The model fitting function to be used. Must take `x` and
#' `y` as first two positional arguments. Further arguments can be passed via
#' `optional_fitter_args` (below). Its return value must be an S3 object with a `plot()` 
#' method.
#' @param save_dir string. The directory in which to store the model and (if `save_plots 
#' == TRUE`) the plots in.
#' @param optional_fitter_args list. Optional arguments passed to `fitter`, e.g. alpha 
#' in case of an elastic net. Default is `list()`, i.e., no arguments other than `x`, `y`
#' passed to `fitter`.
#' @param create_save_dir logical. Whether to create `save_dir` if it does not exist, yet. 
#' Default is `TRUE`.
#' @param save_plots logical. Whether to store a lambda-versus-deviance plot. Default is 
#' `TRUE`.
#' @param plot_fname string. The name of the lambda-versus-deviance plot file inside
#' `save_dir`. Default is `"lambda_vs_deviance.pdf"`.
#' @param fit_fname string. The name of the (binary) model-fit file inside `save_dir`.
#' Default is `"fit_obj.rds"`.
#' @return A `zeroSum.fit` object. The model fit.
#' @details The model fit is stored in 
#' `results_dir`/`save_model_as[1]`/`model`/`save_model_as[2]` with, e.g., `save_model_as 
#' = c("schmitz", "with_ipi")`, meaning that `x` holds gene expression dara and the pheno
#' variabes needed to calculate the IPI score. Row names in `x` and `y` must match.
#' @export
fit <- function(
    x,
    y,
    fitter,
    save_dir,
    optional_fitter_args = list(),
    create_save_dir = TRUE,
    save_plots = TRUE,
    plot_fname = "lambda_vs_deviance.pdf",
    fit_fname = "fit_obj.rds"
){
    check_fitter(fitter, optional_fitter_args)
    # make sure directory to store model exists
    if(!dir.exists(save_dir)){
        if(create_save_dir){
            message("Creating directory ", save_dir, " to store model there")
            dir.create(save_dir, recursive = TRUE)
        } else {
            stop("Path", save_dir, "does not exist and create_path is FALSE")
        }
    }

    fit_obj <- do.call(
        fitter, 
        args = c(list("x" = x, "y" = y), optional_fitter_args)
    )

    # store lambda-vs-deviance plot
    if(save_plots){
        grDevices::pdf(file = file.path(save_dir, plot_fname))
        plot(fit_obj)
        grDevices::dev.off()
    }

    save_fname <- file.path(save_dir, fit_fname)
    saveRDS(fit_obj, file = save_fname)
    message("Done! Model fit stored as ", save_fname)

    return(fit_obj)
}