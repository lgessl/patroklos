ass2d_calculate_2d_metric <- function(self, private, data, model, quiet){

    prep <- model$predict(
        data = data,
        quiet = quiet
    )
    actual <- prep[["actual"]]
    benchmark <- prep[["benchmark"]]
    predicted <- prep[["predicted"]]

    # Prepare for loop
    estimate_list <- list()
    estimate_list[[model$name]] <- predicted
    if(!is.null(self$benchmark))
        estimate_list[[self$benchmark]] <- benchmark
    n_tbls <- length(model$split_index) * length(estimate_list)
    tbl_list <- vector("list", n_tbls)

    j <- 1
    for(i in model$split_index){
        for(estimate_name in names(estimate_list)){
            estimate <- estimate_list[[estimate_name]]
            if(length(table(actual[[i]])) != 2)
                stop("Actual values in split ", i, " are not binary. We have (", 
                    paste(actual[[i]], collapse = ", "), ")")
            if(self$y_metric == "logrank"){
                tbl <- logrank_metric(
                    estimate = estimate[[i]],
                    pheno_tbl = data$pheno_tbl,
                    data = data,
                    y_metric = self$y_metric,
                    x_metric = self$x_metric
                )
            } else if(self$y_metric == "precision_ci"){
                lower_boundary <- estimate_name == model$name
                tbl <- precision_ci(
                    estimate = estimate[[i]],
                    actual = actual[[i]],
                    confidence_level = self$ci_level,
                    y_metric = self$y_metric,
                    x_metric = self$x_metric,
                    lower_boundary = lower_boundary
                )
            } else {
                tbl <- metric_with_rocr(
                    estimate = estimate[[i]],
                    actual = actual[[i]],
                    x_metric = self$x_metric,
                    y_metric = self$y_metric
                )
            }
            names(tbl) <- c(
                self$x_metric, 
                self$y_metric, 
                "cutoff"
            )
            tbl[["split"]] <- i
            tbl[["model"]] <- estimate_name       
            tbl_list[[j]] <- tbl
            j <- j+1
        }
    }
    # Combine all splits to one tibble
    tbl <- dplyr::bind_rows(tbl_list)
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]
    self$data <- tbl
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
        rm_na = c(TRUE, TRUE)
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
        rm_na = c(TRUE, TRUE)
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
    data,
    y_metric = "logrank",
    x_metric = "prevalence"
){
    if(!is.data.frame(pheno_tbl))
        stop("`pheno_tbl` must inherit from `data.frame`")
    if(!"Data" %in% class(data))
        stop("`data` must inherit from `Data`")
    if(!all(names(estimate) %in% pheno_tbl[[data$patient_id_col]]))
        stop("`names(estimate)` must be a subset of `pheno_tbl[[\"data$patient_id_col\"]]`")
    
    # We can no longer tolerate NAs
    estimate <- estimate[!is.na(estimate)]

    # We don't want an empty group
    cutoffs <- estimate[estimate > min(estimate)] |> unique() |> sort()
    prevalence <- numeric(length(cutoffs))
    logrank_p <- numeric(length(cutoffs))
    pheno_mat <- pheno_tbl[, c(data$time_to_event_col, data$event_col)] |> as.matrix()
    rownames(pheno_mat) <- pheno_tbl[[data$patient_id_col]]
    pheno_mat <- pheno_mat[names(estimate), ]
    time <- pheno_mat[, data$time_to_event_col]
    event <- pheno_mat[, data$event_col]
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