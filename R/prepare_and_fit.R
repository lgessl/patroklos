#' @title Prepare data, fit and store models across splits
#' @description Given an expression matrix and a pheno tibble, prepare the data for a 
#' model and fit this model. Do this for all splits into training and test cohort. 
#' Howeve, we do not support multiple time cutoffs at this step; enabling them is 
#' the job of `[training_camp()]`.
#' @param expr_mat numeric matrix. The expression matrix with genes in rows and samples
#' in columns.
#' @param pheno_tbl tibble. The pheno data with samples in rows and variables in columns.
#' @param data_spec DataSpec S3 object. Specifications on the data. See the the 
#' constructor `DataSpec()` for details.
#' @param model_spec ModelSpec S3 object. Specifications on the model. See
#' the the constructor `ModelSpec()` for details.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
#' @return A list of fit objects as returned by the `fit()` method of `model_spec`. 
#' @export
prepare_and_fit <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec,
    quiet = FALSE,
    msg_prefix = ""
){
    if(!inherits(model_spec, "ModelSpec"))
        stop("model_spec_list must be a list of ModelSpec objects")
    if(!inherits(data_spec, "DataSpec"))
        stop("data_spec must be a DataSpec object")

    # Extract
    directory <- model_spec$directory
    # Ensure model directory exists
    if(!dir.exists(directory)){
        message(msg_prefix, "Creating ", directory)
        dir.create(directory, recursive = TRUE)
    }
    if(is.null(data_spec$cohort))
        data_spec$cohort <- "train"
    # Set up list holding fits
    fits <- list()
    stored_fits_fname <- file.path(model_spec$directory, model_spec$fit_fname)
    if(file.exists(stored_fits_fname)){
        if(!quiet) message(msg_prefix, "Found stored fits")
        fits <- readRDS(stored_fits_fname)
    }
    fits[["model_spec"]] <- model_spec
    
    # Fit split after split
    for(i in model_spec$split_index){
        split_name <- paste0(data_spec$split_col_prefix, i)
        if(split_name %in% names(fits)){
            message(msg_prefix, "Found a fit for split ", i, ". Skipping")
            next
        }
        split_ms <- model_spec
        split_ms$split_index <- i
        # Prepare and fit
        x_y <- prepare(
            expr_mat = expr_mat,
            pheno_tbl = pheno_tbl,
            model_spec = split_ms,
            data_spec = data_spec
        )
        x_y <- intersect_by_names(x_y[[1]], x_y[[2]], rm_na = TRUE)
        x <- x_y[[1]]
        y <- x_y[[2]]
        qc_prefit(
            x = x,
            y = y
        )
        fits[[split_name]] <- do.call(
            model_spec$fitter, 
            args = c(
                list("x" = x, "y" = y), 
                model_spec$optional_fitter_args
            )
        )
        message(msg_prefix, "Fitted split ", i, " of ", length(model_spec$split_index))
    }

    saveRDS(fits, stored_fits_fname)

    # Plots about fitting as a grid
    grDevices::pdf(file = file.path(directory, model_spec$plot_fname))
    graphics::par(mfrow = c(ceiling((length(fits)-1)/model_spec$plot_ncols), model_spec$plot_ncols))
    for(i in model_spec$split_index){
        split_name <- paste0(data_spec$split_col_prefix, i)
        fit <- fits[[split_name]]
        plot(fit)
        graphics::title(main = paste0("Split ", i))
    }
    grDevices::dev.off()

    return(fits)
}