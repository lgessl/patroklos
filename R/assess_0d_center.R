assess_0d_center <- function(
    ass_spec_0d,
    data_spec,
    model_spec_list,
    cohorts = c("train", "test"),
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
        res_tbl_list[[cohort]] <- res_tbl
        if(cohort == "test")
            ass_spec_0d$file <- mirror_directory(
                filepath = ass_spec_0d$file,
                mirror = model_tree_mirror
            )
        file <- ass_spec_0d$file
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