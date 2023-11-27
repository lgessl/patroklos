#' @title Fit and store a model
#' @description Given a model and data in a model-specific form, fit this model 
#' and store it
#' @param x numeric matrix. The predictor matrix with samples in rows and their names
#' as row names.
#' @param y numeric matrix. The response matrix with samples in rows and their names
#' as row names.
#' @param model string. The name of the model to fit. Currently supported are
#' * `cox_lasso_zerosum`, Cox proportional hazards with LASSO regularization and
#' zero-sum constraint,
#' * `lasso_zerosum`, ordinary linear regression with LASSO regularization.
#' @param save_model_as character vector of length 2. The first element is the data-set
#' name from which data originates, the second element is the name of the very model
#' (usually a telling name in terms of `x`). See details.
#' @param additional_model_par list. Additional parameters to be passed to the model
#' fitting function, see `zeroSum::zeroSum` for details. Default is `list()`.
#' @param results_dir string. The directory holding results. Default is `"results"`.
#' @param create_dir logical. Whether to create the directory to store the model and 
#' plots in if it does not exist. Default is `TRUE`.
#' @param plots logical. Whether to store a lambda-versus-deviance plot. Default is 
#' `TRUE`.
#' @return A `zeroSum.fit` object. The model fit.
#' @details The model fit is stored in 
#' `results_dir`/`save_model_as[1]`/`model`/`save_model_as[2]` with, e.g., `save_model_as 
#' = c("schmitz", "with_ipi")`, meaning that `x` holds gene expression dara and the pheno
#' variabes needed to calculate the IPI score. Row names in `x` and `y` must match.
#' @export
fit <- function(
    x,
    y,
    model,
    save_model_as,
    additional_model_par = list(),
    results_dir = "results",
    create_dir = TRUE,
    plots = TRUE
){
    # choose model fit function
    fit_func <- NULL
    if(model == "cox_lasso_zerosum"){
        fit_func <- zeroSum::zeroSum
        additional_model_par[["family"]] <- "cox"
    }
    if(model == "lasso_zerosum"){
        fit_func <- zeroSum::zeroSum
        additional_model_par[["family"]] <- "binomial"
    }
    if(is.null(fit_func)){
        stop("Model ", model, " is not supported")
    }

    # make sure directory to store model exists
    save_dir <- file.path(results_dir, save_model_as[1], model, save_model_as[2])
    if(!dir.exists(save_dir)){
        if(create_dir){
            message("Creating directory ", save_dir, " to store model there")
            dir.create(save_dir, recursive = TRUE)
        } else {
            stop("Path", save_dir, "does not exist and create_path is FALSE")
        }
    }
  
    fit_obj <- do.call(
        fit_func, 
        args = c(list("x" = x, "y" = y), additional_model_par)
    )

    # store lambda-vs-deviance plot
    if(plots){
        grDevices::pdf(file = file.path(save_dir, "lambda_vs_deviance.pdf"))
        plot(fit_obj)
        grDevices::dev.off()
    }

    save_fname <- file.path(save_dir, "fit_obj.rds")
    saveRDS(fit_obj, file = save_fname)
    message("Done! Model fit stored as ", save_fname)

    return(fit_obj)
}