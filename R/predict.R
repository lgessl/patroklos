model_predict <- function(self, private, data, quiet){

    # Retrieve model
    model_path <- file.path(self$directory, self$file)
    if(!file.exists(model_path)){
        stop("Model object does not exist at ", model_path)
    }
    fit_obj <- readRDS(model_path)$fit_obj

    benchmark <- NULL
    if(!is.null(data$benchmark_col)){
        benchmark <- data$pheno_tbl[[data$benchmark_col]]
        names(benchmark) <- data$pheno_tbl[[data$patient_id_col]]
    }

    if (data$cohort == "val_predict") {
        y_cox <- as.matrix(data$pheno_tbl[, c(data$time_to_event_col, data$event_col)])
        rownames(y_cox) <- data$pheno_tbl[[data$patient_id_col]]
        yhat <- fit_obj$val_predict
        if (!is.numeric(yhat))
            stop("val_predict attribute of fit object must be a numeric vector")
    } else {
        self$time_cutoffs <- data$pivot_time_cutoff
        max_combo <- 0
        if (!is.null(self$include_from_discrete_pheno))
            max_combo <- fit_obj$combine_n_max_categorical_features
        if (is.null(max_combo))
            stop("No combine_n_max_categorical_features attribute found in ", 
                "fit object")
        self$combine_n_max_categorical_features <- max_combo
        prep <- data$prepare(self, quiet = quiet)
        y_cox <- prep[["y_cox"]]
        yhat <- predict(fit_obj, newx = prep[["x"]])
        rownames(yhat) <- rownames(prep[["x"]])
    }
    y <- binarize_y(y_cox = y_cox, time_cutoff = data$pivot_time_cutoff, 
        pivot_time_cutoff = data$pivot_time_cutoff)
    yhat_y <- intersect_by_names(yhat, y, rm_na = c(TRUE, TRUE))

    res <- list(
        "predicted" = yhat_y[[1]][, 1],
        "actual" = yhat_y[[2]][, 1],
        "benchmark" = benchmark
    )
    return(res)
}