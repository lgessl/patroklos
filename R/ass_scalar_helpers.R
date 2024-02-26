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
    core <- function(i){
        predicted <- prep[["predicted"]][[i]]
        actual <- prep[["actual"]][[i]]
        pa <- intersect_by_names(predicted, actual, rm_na = TRUE)
        do.call(self$metric, list(
            "predicted" =  pa[[1]],
            "actual" = pa[[2]]
            )
        )
    }
    metric <- vapply(seq_along(model$split_index), core, numeric(1))
    return(metric)
}

ass_scalar_assess_center = function(self, private, data, model_list,
    model_tree_mirror, quiet){

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
            "cutoff" = numeric(m),
            "mean" = numeric(m),
            "sd" = numeric(m),
            "min" = numeric(m),
            "max" = numeric(m)
        )
        for(j in 1:m){
            time_cutoff <- model$time_cutoffs[j]
            model_cutoff <- model$at_time_cutoff(time_cutoff)
            res_tbl[j, 1] <- model$name
            res_tbl[j, 2] <- data$cohort
            res_tbl[j, 3] <- time_cutoff
            metric <- self$assess(
                data,
                model = model_cutoff,
                quiet = quiet
            )
            res_tbl[j, 4] <- round(mean(metric), digits = digits) 
            res_tbl[j, 5] <- round(stats::sd(metric), digits = digits)
            res_tbl[j, 6] <- round(min(metric), digits = digits)
            res_tbl[j, 7] <- round(max(metric), digits = digits)
        }
        res_tbl
    }

    res_tbl_list <- list()
    for(cohort in cohorts){
        data$cohort <- cohort
        tbl_list <- lapply(seq_along(model_list), core)
        res_tbl <- dplyr::bind_rows(tbl_list)
        res_tbl <- res_tbl[order(-res_tbl[["mean"]]), ]
        res_tbl_list[[cohort]] <- res_tbl
        file <- self$file
        if(!is.null(file)){
            if(cohort == "test")
                file <- mirror_path(
                    filepath = file,
                    mirror = model_tree_mirror
                )
            if(!dir.exists(dirname(file)))
                dir.create(dirname(file))
            readr::write_csv(res_tbl, file = file)
            if(!quiet)
                message("Writing ", file)
        }
    }
    return(dplyr::bind_rows(res_tbl_list))
}

get_auc <- function(
    predicted,
    actual
){
    pred_obj <- ROCR::prediction(predictions = predicted, labels = actual)
    res <- ROCR::performance(pred_obj, measure = "auc")
    return(res@y.values[[1]])
}