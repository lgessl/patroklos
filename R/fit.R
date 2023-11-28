#' @title Fit and store a model
#' @description Given a model and data in a model-specific form, fit this model 
#' and store it
#' @param x numeric matrix. The predictor matrix with samples in rows and their names
#' as row names.
#' @param y numeric matrix. The response matrix with samples in rows and their names
#' as row names.
#' @param model_spec ModelSpec S3 object. Specifications on the model. See the the 
#' constructor `ModelSpec()` for details.
#' @param save_plots logical. Whether to save the lambda-vs-deviance plot in 
#' `model_spec$save_dir`. Default is `TRUE`.
#' @return A `zeroSum.fit` object, with a `ModelSpec` object appended.
#' @export
fit <- function(
    x,
    y,
    model_spec,
    save_plots = TRUE
){
    # extract values from model_spec
    fitter <- model_spec$fitter
    optional_fitter_args <- model_spec$optional_fitter_args
    save_dir <- model_spec$save_dir
    create_save_dir <- model_spec$create_save_dir
    plot_fname <- model_spec$plot_fname
    fit_fname <- model_spec$fit_fname

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