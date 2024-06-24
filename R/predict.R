model_predict <- function(self, private, data, quiet){

    # Retrieve model
    fit_path <- file.path(self$directory, self$fit_file)
    if(!file.exists(fit_path)){
        stop("Model object does not exist at ", fit_path)
    }
    fits <- readRDS(fit_path)$fits

    predicted_list <- vector("list", length(self$split_index))
    actual_list <- vector("list", length(self$split_index))

    benchmark <- NULL
    benchmark_list <- NULL
    if(!is.null(data$benchmark_col)){
        benchmark <- data$pheno_tbl[[data$benchmark_col]]
        names(benchmark) <- data$pheno_tbl[[data$patient_id_col]]
        benchmark_list <- vector("list", length(self$split_index))
    }

    for(i in self$split_index){
        split_name <- paste0(data$split_col_prefix, i)
        fit <- fits[[split_name]]
        if(is.null(fit))
            stop("No fit found for split ", split)
        if (data$cohort == "val_predict") {
            y_cox <- as.matrix(data$pheno_tbl[, c(data$time_to_event_col, data$event_col)])
            rownames(y_cox) <- data$pheno_tbl[[data$patient_id_col]]
            yhat <- fit$val_predict_list[[fit$lambda_min_index]]
            if (!is.numeric(yhat))
                stop("val_predict attribute of fit object must be a numeric vector")
        } else {
            split_model <- self$clone()
            split_model$split_index <- i
            split_model$time_cutoffs <- data$pivot_time_cutoff
            max_combo <- 0
            if (!is.null(split_model$include_from_discrete_pheno))
                max_combo <- fit$combine_n_max_categorical_features
            if (is.null(max_combo))
                stop("No combine_n_max_categorical_features attribute found in ", 
                    "fit object for split ", split)
            split_model$combine_n_max_categorical_features <- max_combo
            prep <- data$prepare(split_model, quiet = quiet)
            y_cox <- prep[["y_cox"]]
            yhat <- predict(fit, newx = prep[["x"]])
            rownames(yhat) <- rownames(prep[["x"]])
        }
        y <- binarize_y(y_cox = y_cox, time_cutoff = data$pivot_time_cutoff, 
            pivot_time_cutoff = data$pivot_time_cutoff)
        y_yhat <- intersect_by_names(y, yhat, rm_na = c(TRUE, TRUE))
        actual_list[[i]] <- y_yhat[[1]][, 1]
        predicted_list[[i]] <- y_yhat[[2]][, 1]
        if(!is.null(benchmark))
            benchmark_list[[i]] <- benchmark[rownames(y)]
    }

    res <- list(
        "predicted" = predicted_list, 
        "actual" = actual_list,
        "benchmark" = benchmark_list
    )
    return(res)
}