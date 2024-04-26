model_fit <- function(self, private, data, quiet, msg_prefix){

    if(!inherits(data, "Data"))
        stop("data must be a Data object")
    if(is.null(data$cohort)) # Usually fit to train data
        stop("You need to specify the cohort in the Data object")

    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    # Ensure self directory exists
    if(!dir.exists(self$directory)){
        dir.create(self$directory, recursive = TRUE)
        if(!quiet)
            message(msg_prefix, "Creating ", self$directory)
    }
    fits <- list()
    stored_models_file <- file.path(self$directory, self$fit_file)
    if(file.exists(stored_models_file)){
        if(!quiet) message(msg_prefix, "Found stored model")
        fits <- readRDS(stored_models_file)$fits
    }
    
    # Fit split after split
    for(i in self$split_index){
        split_name <- paste0(data$split_col_prefix, i)
        if(split_name %in% names(fits)){
            if(!quiet)
                message(msg_prefix, "Found a fit for split ", i, ". Skipping")
            next
        }
        split_model <- self$clone()
        split_model$split_index <- i
        # Prepare and fit
        x_y <- data$prepare(split_model)
        x_y <- intersect_by_names(x_y[[1]], x_y[[2]], rm_na = TRUE)
        x <- x_y[[1]]
        y <- x_y[[2]]
        qc_prefit(
            x = x,
            y = y
        )
        fits[[split_name]] <- do.call(
            self$fitter, 
            args = c(
                list("x" = x, "y" = y), 
                self$hyperparams,
                list(
                    "n_folds" = self$n_folds, 
                    "append_to_includes" = self$append_to_includes
                )
            )
        )
        if(!quiet)
            message(msg_prefix, "Fitted split ", i, " of ", 
                length(self$split_index), " (", 
                round.POSIXt(Sys.time(), units = "secs"), ")")
    }
    self$fits <- fits

    saveRDS(self, stored_models_file)

    # Plots about fitting as a grid
    grDevices::pdf(file = file.path(self$directory, self$plot_file))
    graphics::par(mfrow = c(1, self$plot_ncols))
    for(i in self$split_index){
        split_name <- paste0(data$split_col_prefix, i)
        fit <- fits[[split_name]]
        plot(fit)
        graphics::title(main = paste0("Split ", i), line = self$plot_title_line)
    }
    grDevices::dev.off()

    invisible(self)
}