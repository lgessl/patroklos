ass_scalar_initialize <- function(
    self, private, metrics, prev_range, benchmark,
    file, round_digits) {
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

ass_scalar_assess <- function(self, private, data, model, quiet) {
    prep <- model$predict(
        data = data,
        quiet = quiet
    )
    if(is.null(prep)) return(rep(NA, length(self$metrics)))
    res_mat <- matrix(numeric(length(self$metrics)), nrow = 1, byrow = TRUE)
    stopifnot(all(names(prep[["predicted"]]) == names(prep[["actual"]])))
    metric_v <- get_metric(
        ass_scalar = self,
        predicted = prep[["predicted"]],
        actual = prep[["actual"]],
        data = data,
        quiet = quiet
    )
    names(metric_v) <- self$metrics
    metric_v
}

ass_scalar_assess_center <- function(
    self, private, data, model_list, mirror, quiet) {
    if (is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()

    # One model corresponds to one row in result tibble
    core <- function(i) {
        model <- model_list[[i]]
        if (!quiet) message("Assessing ", model$name)
        res_tbl <- tibble::tibble(
            "model" = model$name,
            "cohort" = data$cohort
        )
        addon_mat <- matrix(.0, nrow = 1, ncol = length(self$metrics))
        colnames(addon_mat) <- self$metrics
        addon_tbl <- tibble::as_tibble(addon_mat)
        res_tbl <- dplyr::bind_cols(res_tbl, addon_tbl)
        metric_v <- self$assess(
            data,
            model = model,
            quiet = quiet
        )
        res_tbl[1, 3:ncol(res_tbl)] <- as.list(metric_v)
        res_tbl
    }

    tbl_list <- lapply(seq_along(model_list), core)
    res_tbl <- dplyr::bind_rows(tbl_list)
    sortby <- ifelse(length(self$metrics) == 1, "mean", self$metrics[1])
    res_tbl <- res_tbl[order(-res_tbl[[sortby]]), ]
    if (!is.null(self$file)) {
        if (!dir.exists(dirname(self$file))) dir.create(dirname(self$file))
        comment <- paste0(
            "# ", data$name, " | ", data$time_to_event_col, " < ",
            data$pivot_time_cutoff
        )
        header <- paste0(colnames(res_tbl), collapse = ",")
        readr::write_lines(c(comment, header), file = self$file)
        readr::write_csv(res_tbl, file = self$file, append = TRUE)
        if (!quiet) {
            message("Writing ", self$file)
        }
    }
    return(res_tbl)
}

get_metric <- function(
    ass_scalar,
    predicted,
    actual,
    data,
    quiet) {
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
            if (!quiet) {
                message(
                    "No prevalences in the range (",
                    paste(ass_scalar$prev_range, collapse = ", "), "). ",
                    "Available prevalences are (",
                    paste(round(prevs, 3), collapse = ", "), ")."
                )
            }
            return(rep(NA, length(res)))
        }
    }
    swap_sign <- FALSE

    check_one_threshold <- function() {
        if (length(thresholds) != 1) {
            stop(
                ass_scalar$metrics[j], "is only defined for truly binary classfiers ",
                "and it's nonsense to threshold a classfier with ",
                "continous predictions by maximizing the", ass_scalar$metrics[j],
                ". So provide a binary classifier or make sure that ",
                "something reasonable to tune the threshold with (like accuracy) comes ",
                "before \"", ass_scalar$metrics[j], "\" in the `metrics` attribute of the ",
                "AssScalar object."
            )
        }
    }
    for (j in seq_along(ass_scalar$metrics)) {
        switch(ass_scalar$metrics[j],
            "auc" = {
                pred_obj <- ROCR::prediction(predictions = predicted, labels = actual)
                j_res <- ROCR::performance(pred_obj, measure = "auc")@y.values[[1]]
            },
            "accuracy" = {
                j_res <- vapply(
                    thresholds,
                    function(t) mean(as.numeric(predicted >= t) == actual), numeric(1)
                )
            },
            "precision" = {
                j_res <- vapply(
                    thresholds,
                    function(t) mean(actual[predicted >= t]), numeric(1)
                )
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
                j_res <- vapply(
                    thresholds, function(t) {
                        predicted_binary <- predicted >= t
                        if (all(predicted_binary) || all(!predicted_binary)) {
                            return(Inf)
                        }
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
