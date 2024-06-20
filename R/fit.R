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
        prep <- data$prepare(split_model, quiet = quiet)
        qc_prefit(prep)
        fits[[split_name]] <- do.call(
            unitune(self$fitter), 
            args = c(
                list(
                    x = prep[["x"]], 
                    y_bin = prep[["y_bin"]], 
                    y_cox = prep[["y_cox"]], 
                    combine_n_max_categorical_features = self$combine_n_max_categorical_features
                ),
                self$hyperparams
            )
        )
        if(!quiet)
            message(msg_prefix, "Fitted split ", i, " of ", 
                length(self$split_index), " (", 
                round.POSIXt(Sys.time(), units = "secs"), ")")
    }
    self$fits <- fits
    saveRDS(self, stored_models_file)
    invisible(self)
}