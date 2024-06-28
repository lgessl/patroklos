model_fit <- function(self, private, data, quiet, msg_prefix){

    if(!inherits(data, "Data"))
        stop("data must be a Data object")
    if(is.null(data$cohort))
        stop("You need to specify the cohort in the Data object")

    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    if(!dir.exists(self$directory)){
        dir.create(self$directory, recursive = TRUE)
        if(!quiet)
            message(msg_prefix, "Creating ", self$directory)
    }
    model_path <- file.path(self$directory, self$file)
    if(file.exists(model_path)){
        if(!quiet) {
            message(msg_prefix, "Found stored model.")
            if (!is.null(readRDS(model_path)$fit_obj))
                message("Found a fit object in it. Skipping.")
        }
    } else {
        # Prepare and fit
        prep <- data$prepare(model, quiet = quiet)
        qc_prefit(prep)
        model$fit_obj <- do.call(
            unitune(self$fitter), 
            args = c(
                list(
                    x = prep[["x"]], 
                    y_cox = prep[["y_cox"]], 
                    time_cutoffs = model$time_cutoffs,
                    combine_n_max_categorical_features = self$combine_n_max_categorical_features,
                    pivot_time_cutoff = data$pivot_time_cutoff,
                    val_error_fun = self$val_error_fun
                ),
                self$hyperparams
            )
        )
        saveRDS(self, model_path)
    }
    invisible(self)
}