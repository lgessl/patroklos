#' @title Assess multiple models on a single data set with one single value
#' @description Map every model to a 0-dim metric (a single real number) and store 
#' these metrics together with some more analytics across the splits in a sorted 
#' tibble. At the core, we call `[assess_0d()]`.
#' @inheritParams assess_2d_center
#' @param ass_spec_0d AssSpec0d S3 object. See the constructor [`AssSpec0d()`] for more 
#' details.
#' @param model_tree_mirror character vector of length 2. If you want to store the 
#' resulting tibbles, get the test-cohort tibble's file name by mirroring 
#' `ass_spec_0d$file` according to `model_tree_mirror` (see [`mirror_path()`]).
#' @return A tibble. Every row holds the metric (and more analysis across the splits 
#' like standard deviation, minimum, maximum value) for one model in `model_spec_list` 
#' and every time cutoff specified for this model.
#' @export
assess_0d_center <- function(
    ass_spec_0d,
    data_spec,
    model_spec_list,
    cohorts = c("test", "train"),
    model_tree_mirror = c("models", "results"),
    quiet = FALSE
){
    cohorts <- match.arg(cohorts, several.ok = TRUE)
    data <- read(data_spec)
    
    core <- function(i){
        model_spec <- model_spec_list[[i]]
        m <- length(model_spec$time_cutoffs)
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
            time_cutoff <- model_spec$time_cutoffs[j]
            ms_cutoff <- at_time_cutoff(model_spec, time_cutoff)
            res_tbl[j, 1] <- model_spec$name
            res_tbl[j, 2] <- data_spec$cohort
            res_tbl[j, 3] <- time_cutoff
            metric <- assess_0d(
                expr_mat = data[["expr_mat"]],
                pheno_tbl = data[["pheno_tbl"]],
                data_spec = data_spec,
                model_spec = ms_cutoff,
                ass_spec_0d = ass_spec_0d
            )
            res_tbl[j, 4] <- mean(metric)
            res_tbl[j, 5] <- stats::sd(metric)
            res_tbl[j, 6] <- min(metric)
            res_tbl[j, 7] <- max(metric)
        }
        res_tbl
    }

    res_tbl_list <- list()
    for(cohort in cohorts){
        data_spec$cohort <- cohort
        tbl_list <- lapply(seq_along(model_spec_list), core)
        res_tbl <- dplyr::bind_rows(tbl_list)
        res_tbl <- res_tbl[order(res_tbl[["mean"]]), ]
        res_tbl_list[[cohort]] <- res_tbl
        file <- ass_spec_0d$file
        if(cohort == "test")
            file <- mirror_path(
                filepath = file,
                mirror = model_tree_mirror
            )
        if(!is.null(file)){
            if(!dir.exists(dirname(file)))
                dir.create(dirname(file))
            readr::write_csv(res_tbl, file = file)
            if(!quiet)
                message("Writing ", file)
        }
    }
    return(dplyr::bind_rows(res_tbl_list))
}