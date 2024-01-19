calculate_2d_metric <- function(
    actual,
    predicted,
    perf_plot_spec,
    model_spec,
    benchmark = NULL
){
    if(xor(is.null(benchmark), is.null(perf_plot_spec$benchmark))){
        stop("`perf_plot_spec$benchmark` and `benchmark` must be both NULL or both ",
            "not NULL")
    }
    # Extract most frequently used values
    x_metric <- perf_plot_spec$x_metric
    y_metric <- perf_plot_spec$y_metric
    # Prepare for loop
    tbl_list <- list()
    aucs <- rep(NA, length(actual))
    estimate_list <- list()
    estimate_list[[model_spec$name]] <- predicted
    if(!is.null(benchmark)){
        estimate_list[[perf_plot_spec$benchmark]] <- benchmark
    }

    for(i in model_spec$split_index){
        for(estimate_name in names(estimate_list)){
            estimate <- estimate_list[[estimate_name]]
            estimate_actual <- intersect_by_names(
                estimate[[i]], 
                actual[[i]], 
                rm_na = TRUE
            )
            if(length(table(estimate_actual[[2]])) != 2){
                stop("Actual values in split ", i, " are not binary")
            }
            # Calculate performance measures
            rocr_prediction <- ROCR::prediction(
                estimate_actual[[1]], 
                estimate_actual[[2]]
            )
            rocr_perf <- ROCR::performance(
                rocr_prediction,
                measure =  y_metric,
                x.measure = x_metric
            )
            if(estimate_name == model_spec$name){
                # Infer ROC-AUC
                aucs[i] <- ROCR::performance(
                    rocr_prediction,
                    measure = "auc"
                )@y.values[[1]]
            }
            # Store them in a tibble
            tbl <- tibble::tibble(
                rocr_perf@x.values[[1]],
                rocr_perf@y.values[[1]],
                rocr_perf@alpha.values[[1]]
            )
            names(tbl) <- c(
                x_metric, 
                y_metric,
                "cutoff"
            )
            tbl[["split"]] <- i
            tbl[["model"]] <- estimate_name
            tbl_list <- c(tbl_list, list(tbl))
        }
    }
    auc_avg <- mean(aucs, na.rm = TRUE)
    avg_string <- ifelse(length(aucs) > 1, "avg ", "")

    # Combine all splits to one tibble
    tbl <- dplyr::bind_rows(tbl_list)
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]
    perf_plot_spec$data <- tbl
    perf_plot_spec$title <- stringr::str_c(
        perf_plot_spec$title, ", ", 
        "ROC-AUC = ", round(auc_avg, 3)
    )

    return(perf_plot_spec)
}