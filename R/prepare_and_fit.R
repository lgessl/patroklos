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
            print(class(data_spec))
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
            response_type = model_spec$response_type,
            include_from_continuous_pheno = model_spec$include_from_continuous_pheno,
            include_from_discrete_pheno = model_spec$include_from_categorical_pheno,
            patient_id_col = data_spec$patient_id_col,
            pfs_col = data_spec$pfs_col,
            progression_col = data_spec$progression_col,
            pfs_leq = model_spec$pfs_leq
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
    return(fits)
}