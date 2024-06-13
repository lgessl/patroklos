ass_scalar_initialize <- function(self, private, metrics, prev_range, benchmark, 
    file, round_digits){

    stopifnot(is.character(metrics))
    available_metrics <- eval(formals(self$initialize)[["metrics"]])
    stopifnot(all(metrics %in% available_metrics))
    stopifnot(is.numeric(prev_range) && prev_range[1] >= 0 && 
        prev_range[2] <= 1 && prev_range[1] < prev_range[2])
    stopifnot(is.null(benchmark) || is.character(benchmark))
    stopifnot(is.null(file) || is.character(file))
    stopifnot(is.numeric(round_digits) && round_digits >= 0)
    self$metrics <- metrics
    self$prev_range <- prev_range
    self$benchmark <- benchmark
    self$file <- file
    self$round_digits <- round_digits

    invisible(self)
}

ass_scalar_assess <- function(self, private, data, model, quiet){

    prep <- model$predict(
        data = data,
        quiet = quiet
    )
    res_mat <- matrix(.0, nrow = length(model$split_index), 
        ncol = length(self$metrics), byrow = TRUE)
    for (i in seq_along(model$split_index)) {
        predicted <- prep[["predicted"]][[i]]
        actual <- prep[["actual"]][[i]]
        stopifnot(all(names(predicted) == names(actual)))
        res_mat[i, ] <- get_metric(
            ass_scalar = self,
            predicted = predicted,
            actual = actual,
            data = data,
            split_index = i,
            quiet = quiet
        )
    }
    colnames(res_mat) <- self$metrics
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
        if (!quiet) message("Assessing ", model$name)
        m <- length(model$time_cutoffs)
        res_tbl <- tibble::tibble(
            "model" = character(m),
            "cohort" = character(m),
            "cutoff" = numeric(m)
        )
        if (length(self$metrics) == 1) {
            ncol_addon <- 4
            colnames_addon <- c("mean", "sd", "min", "max")
        } else {
            ncol_addon <- length(self$metrics)
            colnames_addon <- self$metrics
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
            if (length(self$metrics) == 1) {
                metric <- metric_mat[, 1]
                metric <- metric[!is.na(metric)]
                if (length(metric) == 0) {
                    res_tbl[j, 4:7] <- NA
                    next
                }
                res_tbl[j, 4] <- round(mean(metric), digits = digits) 
                res_tbl[j, 5] <- round(stats::sd(metric), digits = digits)
                res_tbl[j, 6] <- round(min(metric), digits = digits)
                res_tbl[j, 7] <- round(max(metric), digits = digits)
            } else {
                for(k in seq_along(self$metrics)){
                    metric <- metric_mat[, k]
                    metric <- metric[!is.na(metric)]
                    if (length(metric) == 0) {
                        res_tbl[j, 3 + k] <- NA
                        next
                    }
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
        sortby <- ifelse(length(self$metrics) == 1, "mean", self$metrics[1])
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
            comment <- paste0("# ", data$name, ", ", data$time_to_event_col, " < ", 
                data$pivot_time_cutoff)
            header <- paste0(colnames(res_tbl), collapse = ",")
            readr::write_lines(c(comment, header), file = file)
            readr::write_csv(res_tbl, file = file, append = TRUE)
            if(!quiet)
                message("Writing ", file)
        }
    }
    data$cohort <- cohorts # No side effects
    return(dplyr::bind_rows(res_tbl_list))
}

get_metric <- function(
    ass_scalar,
    predicted,
    actual,
    data,
    split_index,
    quiet
){
    # Some metrics want binary data, i.e. in {0, 1} 
    binary_bool <- all(predicted %in% c(0, 1))
    res <- numeric(length(ass_scalar$metrics))
    thresholds <- 1
    if (!binary_bool) {
        thresholds <- unique(predicted)
        prevs <- vapply(thresholds, function(t) mean(predicted >= t), numeric(1))
        thresholds <- thresholds[prevs >= ass_scalar$prev_range[1] & 
            prevs <= ass_scalar$prev_range[2]]
        if (length(thresholds) == 0) {
            if (!quiet)
                message("In split ", split_index, ": No prevalences in the range (", 
                    paste(ass_scalar$prev_range, collapse = ", "), "). ",
                    "Available prevalences are (", 
                    paste(round(prevs, 3), collapse = ", "), ").")
            return(rep(NA, length(res)))
        }
    }
    swap_sign <- FALSE

    check_one_threshold <- function(){
        if(length(thresholds) != 1)
            stop(ass_scalar$metrics[j], "is only defined for truly binary classfiers ",
                "and it's nonsense to threshold a classfier with ", 
                "continous predictions by maximizing the", ass_scalar$metrics[j], 
                ". So provide a binary classifier or make sure that ",
                "something reasonable to tune the threshold with (like accuracy) comes ",
                "before \"", ass_scalar$metrics[j], "\" in the `metrics` attribute of the ",
                "AssScalar object.")
    }
    for (j in seq_along(ass_scalar$metrics)) {
        switch(ass_scalar$metrics[j],
            "auc" = {
                pred_obj <- ROCR::prediction(predictions = predicted, labels = actual)
                j_res <- ROCR::performance(pred_obj, measure = "auc")@y.values[[1]]
            },
            "accuracy" = {
                j_res <- vapply(thresholds, 
                    function(t) mean(as.numeric(predicted >= t) == actual), numeric(1))
            },
            "precision" = {
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
                check_one_threshold()
                j_res <- thresholds
            },
            stop("Unknown metric: ", ass_scalar$metrics[j])
        )
        max_idx <- which.max(ifelse(swap_sign, -1, 1) * j_res)
        swap_sign <- FALSE
        if (length(max_idx) == 1) {
            thresholds <- thresholds[max_idx]
            res[j] <- j_res[max_idx]
        } else {
            res[j] <- NA
        }
    }
    res
}

#' @title Prepend a prefix to the file attribute of a list of AssScalar objects
#' @description Modify the AssScalar objects *in place*. This function supports 
#' the paradigm of specifying the `Model`s for 
#' for *all* data sets at one place and them quickly adopting for a specific data 
#' set.
#' @param ass_scalar_list list of AssScalar objects.
#' @param prefix character Prepend this prefix to the `file` attribute of each 
#' element in `ass_scalar_list`.
#' @return list of modified `AssScalar`s.
#' @export
prepend_to_filename <- function(ass_scalar_list, prefix) {
    stopifnot(is.character(prefix))
    lapply(ass_scalar_list, function(as) {
        stopifnot(inherits(as, "AssScalar"))
        as$file <- file.path(prefix, as$file)
    }) 
    invisible(NULL)
}