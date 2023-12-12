#' @title Prepare data, fit and store models
#' @description Given an expression matrix and a pheno tibble, prepare the data for a 
#' list of models and fit these models to the data.
#' @param expr_mat numeric matrix. The expression matrix with genes in rows and samples
#' in columns.
#' @param pheno_tbl tibble. The pheno data with samples in rows and variables in columns.
#' @param data_spec DataSpec S3 object. Specifications on the data. See the the 
#' constructor `DataSpec()` for details.
#' @param model_spec_list list of ModelSpec S3 objects. Specifications on the models. See
#' the the constructor `ModelSpec()` for details.
#' @return A list of `zeroSum.fit` objects, with a `ModelSpec` object appended.
#' @export
prepare_and_fit <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec_list
){
    fits <- list()
    for(i in seq_along(model_spec_list)){
        model_spec <- model_spec_list[[i]]

        # Checks
        if(!inherits(model_spec, "ModelSpec")){
            stop("model_spec_list must be a list of ModelSpec objects")
        }
        if(!inherits(data_spec, "DataSpec")){
            stop("data_spec must be a DataSpec object")
        }

        model_obj_file <- file.path(model_spec$save_dir, model_spec$fit_fname)
        if(file.exists(model_obj_file)){
            message("Model object already exists at ", model_obj_file, ". Skipping.")
            fits[[i]] <- readRDS(model_obj_file)
            next
        }
        # Prepare in model-specific manner ...
        x_y <- prepare(
            expr_mat = expr_mat,
            pheno_tbl = pheno_tbl,
            model_spec = model_spec,
            data_spec = data_spec
        )
        # ... and fit
        fit_obj <- fit(
            x = x_y[["x"]],
            y = x_y[["y"]],
            model_spec = model_spec
        )
        fit_obj$model_spec <- model_spec
        fits[[i]] <- fit_obj
    }
    invisible(fits)
}