#' @title Prepare data, fit and store models
#' @description Given an expression matrix and a pheno tibble, prepare the data for a 
#' list of models and fit these models to the data.
#' @param expr_mat numeric matrix. The expression matrix with genes in rows and samples
#' in columns.
#' @param pheno_tbl tibble. The pheno data with samples in rows and variables in columns.
#' @param data_spec DataSpec S3 object. Specifications on the data. See the the 
#' constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model. See
#' the the constructor `ModelSpec()` for details.
#' @return A list of fit objects as returned by the `fit()` method of `model_spec`. 
#' @export
prepare_and_fit <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec
){
    if(!inherits(model_spec, "ModelSpec"))
        stop("model_spec_list must be a list of ModelSpec objects")
    if(!inherits(data_spec, "DataSpec"))
        stop("data_spec must be a DataSpec object")

    message(model_spec$name, " on ", data_spec$name, ": begin fitting")
    # Extract
    directory <- model_spec$directory
    # Ensure model directory exists
    if(!dir.exists(directory)){
        message("\tCreating ", directory)
        dir.create(directory, recursive = TRUE)
    }
    # Set up list holding fits
    fits <- list()
    stored_fits_fname <- file.path(model_spec$directory, model_spec$fit_fname)
    if(file.exists(stored_fits_fname)){
        message("Found stored fits")
        fits <- readRDS(stored_fits_fname)
    }
    fits[["model_spec"]] <- model_spec
    
    # Fit split after split
    for(i in model_spec$split_index){
        split_colname <- paste0(data_spec$split_col_prefix, i)
        if(split_colname %in% names(fits)){
            message("\tFound a fit for split ", i, ". Skipping")
            next
        }
        # Reduce data to train samples for split i
        if(!split_colname %in% colnames(pheno_tbl))
            stop("Column ", split_colname, " not found in pheno table.")
        train_bool <- pheno_tbl[[split_colname]] == "train"
        pheno_tbl_train <- pheno_tbl[train_bool, ]
        expr_mat_train <- expr_mat[train_bool, ]
        # Prepare and fit
        x_y <- prepare(
            expr_mat = expr_mat_train,
            pheno_tbl = pheno_tbl_train,
            model_spec = model_spec,
            data_spec = data_spec
        )
        fit_obj <- do.call(
            model_spec$fitter, 
            args = c(
                list("x" = x_y[["x"]], "y" = x_y[["y"]]), 
                model_spec$optional_fitter_args
            )
        )
        fits[[split_colname]] <- fit_obj
    }

    saveRDS(fits, stored_fits_fname)

    # Plots about fitting as a grid
    grDevices::jpeg(file = file.path(directory, model_spec$plot_fname))
    par(mfrow = c(ceiling(length(fits)/model_spec$plot_ncols), model_spec$plot_ncols))
    for(fit in fits)
        if("ModelSpec" %in% class(fit)) next
        plot(fit)
    grDevices::dev.off()

    return(fits)
}