model_predict <- function(self, private, data, quiet){

    if(length(self$time_cutoffs) > 1L)
        stop("Multiple time cutoffs are not supported")
    if(is.null(data$cohort))
        stop("You need to specify the cohort in the Data object")

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
        split_model <- self$clone()
        split_model$split_index <- i
        split_model$time_cutoffs <- data$pivot_time_cutoff
        prep <- data$prepare(split_model, quiet = quiet)
        actual <- prep[["y_bin"]][, 1]
        fit <- fits[[split_name]]
        if(is.null(fit))
            stop("No fit found for split ", split)
        predicted <- predict(fit, newx = prep[["x"]])

        # Check what predict method did
        if(!is.numeric(predicted)){
            stop("predict method for class ", class(fit), " does not return a ", 
            "numeric matrix or vector. Instead it has class ", class(predicted), 
            ".")
        }
        if(is.matrix(predicted)){
            if(ncol(predicted) > 1L){
                stop("predict method for class ", class(fit), " returns a matrix ",
                "with more than one column")
            }
            predicted <- predicted[, 1]
        }
        if(is.null(names(predicted))){
            names(predicted) <- rownames(x_y[["x"]])
        }
        
        predicted_list[[i]] <- predicted
        actual_list[[i]] <- actual
        if(!is.null(benchmark))
            benchmark_list[[i]] <- benchmark[names(actual)]
    }

    res <- list(
        "predicted" = predicted_list, 
        "actual" = actual_list,
        "benchmark" = benchmark_list
    )
    return(res)
}