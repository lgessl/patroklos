calculate_2d_metric <- function(
    actual,
    predicted,
    perf_plot_spec,
    model_spec,
    benchmark = NULL,
    pheno_tbl = NULL,
    data_spec = NULL
){
    if(xor(is.null(benchmark), is.null(perf_plot_spec$benchmark))){
        stop("`perf_plot_spec$benchmark` and `benchmark` must be both NULL or both ",
            "not NULL")
    }
    # Prepare for loop
    tbl_list <- list()
    estimate_list <- list()
    estimate_list[[model_spec$name]] <- predicted
    if(!is.null(benchmark)){
        estimate_list[[perf_plot_spec$benchmark]] <- benchmark
    }

    for(i in model_spec$split_index){
        for(estimate_name in names(estimate_list)){
            estimate <- estimate_list[[estimate_name]]
            if(length(table(actual[[i]])) != 2)
                stop("Actual values in split ", i, " are not binary")
            if(perf_plot_spec$y_metric == "logrank"){
                tbl <- logrank_metric(
                    estimate = estimate[[i]],
                    pheno_tbl = pheno_tbl,
                    data_spec = data_spec,
                    y_metric = perf_plot_spec$y_metric,
                    x_metric = perf_plot_spec$x_metric
                )
            } else {
                tbl <- metric_with_rocr(
                    estimate = estimate[[i]],
                    actual = actual[[i]],
                    x_metric = perf_plot_spec$x_metric,
                    y_metric = perf_plot_spec$y_metric
                )
            }
            tbl[["split"]] <- i
            tbl[["model"]] <- estimate_name       
            tbl_list <- c(tbl_list, list(tbl))
        }
    }
    # Combine all splits to one tibble
    tbl <- dplyr::bind_rows(tbl_list)
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]
    perf_plot_spec$data <- tbl

    return(perf_plot_spec)
}


metric_with_rocr <- function(
    estimate,
    actual,
    x_metric,
    y_metric
){
    estimate_actual <- intersect_by_names(
        estimate, 
        actual, 
        rm_na = TRUE
    )
    # Calculate performance measures
    rocr_prediction <- ROCR::prediction(
        predictions = estimate_actual[[1]], 
        labels = estimate_actual[[2]]
    )
    rocr_perf <- ROCR::performance(
        rocr_prediction,
        measure =  y_metric,
        x.measure = x_metric
    )
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
    return(tbl)
}


logrank_metric <- function(
    estimate,
    pheno_tbl,
    data_spec,
    y_metric = "logrank",
    x_metric = "prevalence"
){
    if(!is.data.frame(pheno_tbl))
        stop("`pheno_tbl` must inherit from `data.frame`")
    if(!"DataSpec" %in% class(data_spec))
        stop("`data_spec` must inherit from `DataSpec`")
    if(!all(names(estimate) %in% pheno_tbl[[data_spec$patient_id_col]]))
        stop("`names(estimate)` must be a subset of `pheno_tbl[[\"data_spec$patient_id_col\"]]`")
    
    # We can no longer tolerate NAs
    estimate <- estimate[!is.na(estimate)]

    # We don't want an empty group
    cutoffs <- estimate[estimate > min(estimate)] |> unique() |> sort()
    prevalence <- numeric(length(cutoffs))
    logrank_p <- numeric(length(cutoffs))
    pheno_mat <- pheno_tbl[, c(data_spec$time_to_event_col, data_spec$event_col)] |> as.matrix()
    rownames(pheno_mat) <- pheno_tbl[[data_spec$patient_id_col]]
    pheno_mat <- pheno_mat[names(estimate), ]
    time <- pheno_mat[, data_spec$time_to_event_col]
    event <- pheno_mat[, data_spec$event_col]
    if(!identical(names(time), names(event)) || !identical(names(time), names(estimate)))
        stop("Names of `time`, `event`, and `estimate` must be identical")

    for(i in seq_along(cutoffs)){
        cutoff <- cutoffs[i]
        groups <- ifelse(estimate >= cutoff, 1, 0)
        prevalence[i] <- mean(groups)
        res <- survival::survdiff(
            survival::Surv(time = time, event = event) ~ groups
        )
        logrank_p[i] <- res[["pvalue"]]
    }
    tbl <- tibble::tibble(
        prevalence,
        logrank_p,
        cutoffs
    )
    names(tbl) <- c(
        x_metric, 
        y_metric,
        "cutoff"
    )
    return(tbl)
}