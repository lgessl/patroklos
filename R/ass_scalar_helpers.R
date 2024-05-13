ass_scalar_initialize <- function(self, private, metric, pivot_time_cutoff, 
    lambda, benchmark, file, round_digits){

    stopifnot(is.character(metric))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(pivot_time_cutoff > 0)
    stopifnot(is.character(lambda) || is.numeric(lambda))
    stopifnot(is.null(benchmark) || is.character(benchmark))
    stopifnot(is.null(file) || is.character(file))
    stopifnot(is.numeric(round_digits) && round_digits >= 0)
    self$metric <- metric
    self$pivot_time_cutoff <- pivot_time_cutoff
    self$lambda <- lambda
    self$benchmark <- benchmark
    self$file <- file
    self$round_digits <- round_digits

    invisible(self)
}

ass_scalar_assess <- function(self, private, data, model, quiet){

    prep <- model$predict(
        data = data, 
        lambda = self$lambda,
        pivot_time_cutoff = self$pivot_time_cutoff
    )
    res_mat <- matrix(.0, nrow = length(model$split_index), 
        ncol = length(self$metric), byrow = TRUE)
    for (i in seq_along(model$split_index)) {
        predicted <- prep[["predicted"]][[i]]
        actual <- prep[["actual"]][[i]]
        pa <- intersect_by_names(predicted, actual, rm_na = TRUE)
        res_mat[i, ] <- get_metric(
            metrics = self$metric,
            predicted = pa[[1]],
            actual = pa[[2]],
            data = data
        )
    }
    colnames(res_mat) <- self$metric
    rownames(res_mat) <- self$split_index
    res_mat
}

ass_scalar_assess_center <- function(self, private, data, model_list,
    mirror, quiet){

    cohorts <- data$cohort
    if(is.null(cohorts)) cohorts <- "test"
    digits <- self$round_digits

    data$read()
    
    core <- function(i){
        model <- model_list[[i]]
        m <- length(model$time_cutoffs)
        res_tbl <- tibble::tibble(
            "model" = character(m),
            "cohort" = character(m),
            "cutoff" = numeric(m)
        )
        if (length(self$metric) == 1) {
            ncol_addon <- 4
            colnames_addon <- c("mean", "sd", "min", "max")
        } else {
            ncol_addon <- length(self$metric)
            colnames_addon <- self$metric
        }
        addon_mat <- matrix(.0, nrow = m, ncol = ncol_addon)
        colnames(addon_mat) <- colnames_addon
        addon_tbl <- tibble::as_tibble(addon_mat)
        res_tbl <- dplyr::bind_cols(res_tbl, addon_tbl)
        for(j in 1:m){
            time_cutoff <- model$time_cutoffs[j]
            model_cutoff <- model$at_time_cutoff(time_cutoff)
            res_tbl[j, 1] <- model$name
            res_tbl[j, 2] <- data$cohort
            res_tbl[j, 3] <- time_cutoff
            metric_mat <- self$assess(
                data,
                model = model_cutoff,
                quiet = quiet
            )
            if (length(self$metric) == 1) {
                metric <- metric_mat[, 1]
                res_tbl[j, 4] <- round(mean(metric), digits = digits) 
                res_tbl[j, 5] <- round(stats::sd(metric), digits = digits)
                res_tbl[j, 6] <- round(min(metric), digits = digits)
                res_tbl[j, 7] <- round(max(metric), digits = digits)
            } else {
                for(k in seq_along(self$metric)){
                    metric <- metric_mat[, k]
                    res_tbl[j, 3 + k] <- round(mean(metric), digits = digits)
                }
            }
        }
        res_tbl
    }

    res_tbl_list <- list()
    for(cohort in cohorts){
        data$cohort <- cohort
        tbl_list <- lapply(seq_along(model_list), core)
        res_tbl <- dplyr::bind_rows(tbl_list)
        sortby <- ifelse(length(self$metric) == 1, "mean", self$metric[1])
        res_tbl <- res_tbl[order(-res_tbl[[sortby]]), ]
        res_tbl_list[[cohort]] <- res_tbl
        file <- self$file
        if(!is.null(file)){
            if(cohort == "test")
                file <- mirror_path(
                    filepath = file,
                    mirror = mirror
                )
            if(!dir.exists(dirname(file)))
                dir.create(dirname(file))
            readr::write_csv(res_tbl, file = file)
            if(!quiet)
                message("Writing ", file)
        }
    }
    data$cohort <- cohorts # No side effects
    return(dplyr::bind_rows(res_tbl_list))
}

get_metric <- function(
    metrics,
    predicted,
    actual,
    data
){
    # Some metrics want binary data, i.e. in {0, 1} 
    binary_bool <- all(predicted %in% c(0, 1))
    res <- numeric(length(metrics))
    thresholds <- ifelse(binary_bool, 1, unique(predicted))
    swap_sign <- FALSE

    check_one_threshold <- function(){
        if(length(thresholds) != 1)
            stop(metrics[j], "is only defined for truly binary classfiers",
                " and it's nonsense to threshold a classfier with ", 
                " continous predictions by maximizing the", metrics[j], 
                ". So provide a binary classifier or make sure that ",
                "something reasonable to tune the threshold with (like accuracy) comes ",
                "before \"", metrics[j], "\" in the `metric` attribute of the ",
                "AssScalar object.")
    }
    for (j in seq_along(metrics)) {
        switch(metrics[j],
            "auc" = {
                pred_obj <- ROCR::prediction(predictions = predicted, labels = actual)
                j_res <- ROCR::performance(pred_obj, measure = "auc")@y.values[[1]]
            },
            "accuracy" = {
                j_res <- vapply(thresholds, 
                    function(t) mean(as.numeric(predicted >= t) == actual), numeric(1))
            },
            "precision" = {
                check_one_threshold()
                j_res <- vapply(thresholds, 
                    function(t) mean(actual[predicted >= t]), numeric(1))
            },
            "prevalence" = {
                check_one_threshold()
                j_res <- mean(predicted >= thresholds)
            },
            "n_true" = {
                check_one_threshold()
                j_res <- sum(actual)
            },
            "perc_true" = {
                check_one_threshold()
                j_res <- mean(actual)
            },
            "n_samples" = {
                check_one_threshold()
                j_res <- length(actual)
            },
            "logrank" = {
                pheno_tbl <- data$pheno_tbl
                pheno_tbl <- pheno_tbl[pheno_tbl[[data$patient_id_col]] %in% names(predicted), ]
                time <- pheno_tbl[[data$time_to_event_col]]
                event <- pheno_tbl[[data$event_col]]
                j_res <- vapply(thresholds, function(t) {
                    predicted_binary <- predicted >= t
                    if (all(predicted_binary) || all(!predicted_binary))
                        return(Inf)
                    survival::survdiff(
                        survival::Surv(time = time, event = event) ~ predicted_binary
                    )$pvalue
                    }, 
                    numeric(1)
                )
                swap_sign <- TRUE
            },
            "threshold" = {
                if (length(thresholds) != 1)
                    stop("Either assess a binary classifier or make sure", 
                        " something reasonable to tune the threshold with",
                        " (like accuracy) comes before \"threshold\" in the", 
                        " `metric` attribute of the AssScalar object.")
                j_res <- thresholds
            },
            stop("Unknown metric: ", metrics[j])
        )
        max_idx <- which.max(ifelse(swap_sign, -1, 1) * j_res)
        thresholds <- thresholds[max_idx]
        res[j] <- j_res[max_idx]
        swap_sign <- FALSE
    }
    res
}
