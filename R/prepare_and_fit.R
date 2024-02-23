#' @title Prepare data, fit and store models across splits
#' @description Given an expression matrix and a pheno tibble, prepare the data for a 
#' model and fit this model. Do this for all splits into training and test cohort. 
#' Howeve, we do not support multiple time cutoffs at this step; enabling them is 
#' the job of `[training_camp()]`.
#' @param expr_mat numeric matrix. The expression matrix with genes in rows and samples
#' in columns.
#' @param pheno_tbl tibble. The pheno data with samples in rows and variables in columns.
#' @param data DataSpec S3 object. Specifications on the data. See the the 
#' constructor `DataSpec()` for details.
#' @param model ModelSpec S3 object. Specifications on the model. See
#' the the constructor `ModelSpec()` for details.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
#' @return A list of fit objects as returned by the `fit()` method of `model`. 
#' @importFrom graphics plot
#' @export
prepare_and_fit <- function(
    expr_mat,
    pheno_tbl,
    data,
    model,
    quiet = FALSE,
    msg_prefix = ""
){
    if(!inherits(model, "ModelSpec"))
        stop("model_list must be a list of ModelSpec objects")
    if(!inherits(data, "DataSpec"))
        stop("data must be a DataSpec object")
    if(is.null(data$cohort)) # Usually fit to train data
        data$cohort <- "train"
    # Extract
    directory <- model$directory
    # Ensure model directory exists
    if(!dir.exists(directory)){
        dir.create(directory, recursive = TRUE)
        if(!quiet)
            message(msg_prefix, "Creating ", directory)
    }
    if(is.null(data$cohort))
        data$cohort <- "train"
    # Set up list holding fits
    fits <- list()
    stored_fits_file <- file.path(model$directory, model$fit_file)
    if(file.exists(stored_fits_file)){
        if(!quiet) message(msg_prefix, "Found stored fits")
        fits <- readRDS(stored_fits_file)
    }
    fits[["model"]] <- model
    
    # Fit split after split
    for(i in model$split_index){
        split_name <- paste0(data$split_col_prefix, i)
        if(split_name %in% names(fits)){
            if(!quiet)
                message(msg_prefix, "Found a fit for split ", i, ". Skipping")
            next
        }
        split_ms <- model
        split_ms$split_index <- i
        # Prepare and fit
        x_y <- prepare(
            expr_mat = expr_mat,
            pheno_tbl = pheno_tbl,
            model = split_ms,
            data = data
        )
        x_y <- intersect_by_names(x_y[[1]], x_y[[2]], rm_na = TRUE)
        x <- x_y[[1]]
        y <- x_y[[2]]
        qc_prefit(
            x = x,
            y = y
        )
        fits[[split_name]] <- do.call(
            model$fitter, 
            args = c(
                list("x" = x, "y" = y), 
                model$optional_fitter_args
            )
        )
        if(!quiet)
            message(msg_prefix, "Fitted split ", i, " of ", 
                length(model$split_index), " (", 
                round.POSIXt(Sys.time(), units = "secs"), ")")
    }

    saveRDS(fits, stored_fits_file)

    # Plots about fitting as a grid
    grDevices::pdf(file = file.path(directory, model$plot_file))
    graphics::par(mfrow = c(1, model$plot_ncol))
    for(i in model$split_index){
        split_name <- paste0(data$split_col_prefix, i)
        fit <- fits[[split_name]]
        plot(fit)
        graphics::title(main = paste0("Split ", i), line = model$plot_title_line)
    }
    grDevices::dev.off()

    return(fits)
}