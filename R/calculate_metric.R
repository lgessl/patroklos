calculate_2d_metric <- function(
    actual,
    predicted,
    ass_spec_2d,
    model_spec,
    benchmark = NULL,
    pheno_tbl = NULL,
    data_spec = NULL
){
    if(xor(is.null(benchmark), is.null(ass_spec_2d$benchmark))){
        stop("`ass_spec_2d$benchmark` and `benchmark` must be both NULL or both ",
            "not NULL")
    }
    # Prepare for loop
    tbl_list <- list()
    estimate_list <- list()
    estimate_list[[model_spec$name]] <- predicted
    if(!is.null(benchmark)){
        estimate_list[[ass_spec_2d$benchmark]] <- benchmark
    }

    for(i in model_spec$split_index){
        for(estimate_name in names(estimate_list)){
            estimate <- estimate_list[[estimate_name]]
            if(length(table(actual[[i]])) != 2)
                stop("Actual values in split ", i, " are not binary")
            if(ass_spec_2d$y_metric == "logrank"){
                tbl <- logrank_metric(
                    estimate = estimate[[i]],
                    pheno_tbl = pheno_tbl,
                    data_spec = data_spec,
                    y_metric = ass_spec_2d$y_metric,
                    x_metric = ass_spec_2d$x_metric
                )
            } else if(ass_spec_2d$y_metric == "precision_ci"){
                lower_boundary <- estimate_name == model_spec$name
                tbl <- precision_ci(
                    estimate = estimate[[i]],
                    actual = actual[[i]],
                    confidence_level = ass_spec_2d$ci_level,
                    y_metric = ass_spec_2d$y_metric,
                    x_metric = ass_spec_2d$x_metric,
                    lower_boundary = lower_boundary
                )
            } else {
                tbl <- metric_with_rocr(
                    estimate = estimate[[i]],
                    actual = actual[[i]],
                    x_metric = ass_spec_2d$x_metric,
                    y_metric = ass_spec_2d$y_metric
                )
            }
            names(tbl) <- c(
                ass_spec_2d$x_metric, 
                ass_spec_2d$y_metric, 
                "cutoff"
            )
            tbl[["split"]] <- i
            tbl[["model"]] <- estimate_name       
            tbl_list <- c(tbl_list, list(tbl))
        }
    }
    # Combine all splits to one tibble
    tbl <- dplyr::bind_rows(tbl_list)
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]
    ass_spec_2d$data <- tbl

    return(ass_spec_2d)
}


precision_ci <- function(
    estimate,
    actual,
    confidence_level = 0.95,
    y_metric = "ci_boundary",
    x_metric = "prevalence",
    lower_boundary = TRUE
){
    estimate_actual <- intersect_by_names(
        estimate, 
        actual, 
        rm_na = TRUE
    )
    estimate <- estimate_actual[[1]]
    actual <- estimate_actual[[2]]
    cutoffs <- estimate[estimate > min(estimate)] |> unique() |> sort()

    precision_ci_core <- function(
        cutoff
    ){
        positive <- ifelse(estimate >= cutoff, 1, 0)
        prevalence <- mean(positive)
        htest <- stats::binom.test(
            x = sum(positive * actual),
            n = sum(positive),
            conf.level = confidence_level
        )
        ci_boundary <- ifelse(lower_boundary, htest$conf.int[1], htest$conf.int[2])
        c("prevalence" = prevalence, "ci_boundary" = ci_boundary)
    }
    mat <- vapply(cutoffs, precision_ci_core, numeric(2))
    if(ncol(mat) == 0){
        mat <- matrix(nrow = 2, ncol = 0)
        rownames(mat) <- c("prevalence", "ci_boundary")
    }

    # Store them in a tibble
    tbl <- tibble::tibble(
        mat["prevalence", ],
        mat["ci_boundary", ],
        cutoffs
    )
    return(tbl)
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
    any_row_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_row_na, ]
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

    logrank_core <- function(
        cutoff
    ){
        groups <- ifelse(estimate >= cutoff, 1, 0)
        prevalence <- mean(groups)
        pval <- survival::survdiff(
            survival::Surv(time = time, event = event) ~ groups
        )[["pvalue"]]
        c("prevalence" = prevalence, "pval" = pval)
    }
    mat <- vapply(cutoffs, logrank_core, numeric(2))
    if(ncol(mat) == 0){
        mat <- matrix(nrow = 2, ncol = 0)
        rownames(mat) <- c("prevalence", "pval")
    }
    tbl <- tibble::tibble(
        mat["prevalence", ],
        mat["pval", ],
        cutoffs
    )
    return(tbl)
}