ass2d_get_2d_metric <- function(self, private, data, model, quiet){

    prep <- model$predict(
        data = data,
        quiet = quiet
    )
    if (is.null(prep)) return(NULL)
    actual <- prep[["actual"]]
    predicted <- prep[["predicted"]]

    if( self$y_metric == "logrank" ){
        tbl <- logrank_metric(
            predicted = predicted,
            pheno_tbl = data$pheno_tbl,
            data = data
        )
    } else if( self$y_metric == "precision_ci" ){
        tbl <- precision_ci(
            predicted = predicted,
            actual = actual,
            confidence_level = self$ci_level
        )
    } else if (self$x_metric == "rank" && self$y_metric == "risk score") {
        tbl <- scores_rank(
           predicted = predicted,
           actual = actual
        )
    } else {
        tbl <- metric_with_rocr(
            predicted = predicted,
            actual = actual,
            x_metric = self$x_metric,
            y_metric = self$y_metric
        )
    }
    names(tbl) <- c(
        self$x_metric, 
        self$y_metric, 
        "cutoff"
    )
    tbl[["model"]] <- model$name
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]
    # Constraint data to range
    tbl <- tbl[
        tbl[[self$x_metric]] >= self$xlim[1] &
        tbl[[self$x_metric]] <= self$xlim[2],
    ]
    tbl <- tbl[
        tbl[[self$y_metric]] >= self$ylim[1] &
        tbl[[self$y_metric]] <= self$ylim[2],
    ]
    return(tbl)
}

precision_ci <- function(
    predicted,
    actual,
    confidence_level
){
    cutoffs <- predicted[predicted > min(predicted)] |> unique() |> sort()

    precision_ci_core <- function(
        cutoff
    ){
        positive <- 1 * (predicted >= cutoff)
        prevalence <- mean(positive)
        htest <- stats::binom.test(
            x = sum(positive * actual),
            n = sum(positive),
            conf.level = confidence_level
        )
        c("prevalence" = prevalence, "ci_boundary" = htest$conf.int[1])
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
    predicted,
    actual,
    x_metric,
    y_metric
){
    # Calculate performance measures
    rocr_prediction <- ROCR::prediction(
        predictions = predicted, 
        labels = actual
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
    predicted,
    pheno_tbl,
    data
){
    if(!is.data.frame(pheno_tbl))
        stop("`pheno_tbl` must inherit from `data.frame`")
    if(!"Data" %in% class(data))
        stop("`data` must inherit from `Data`")
    if(!all(names(predicted) %in% pheno_tbl[[data$patient_id_col]]))
        stop("`names(predicted)` must be a subset of `pheno_tbl[[\"data$patient_id_col\"]]`")

    # We don't want an empty group
    cutoffs <- predicted[predicted > min(predicted)] |> unique() |> sort()
    prevalence <- numeric(length(cutoffs))
    logrank_p <- numeric(length(cutoffs))
    pheno_mat <- pheno_tbl[, c(data$time_to_event_col, data$event_col)] |> as.matrix()
    rownames(pheno_mat) <- pheno_tbl[[data$patient_id_col]]
    pheno_mat <- pheno_mat[names(predicted), ]
    time <- pheno_mat[, data$time_to_event_col]
    event <- pheno_mat[, data$event_col]
    if(!identical(names(time), names(event)) || !identical(names(time), names(predicted)))
        stop("Names of `time`, `event`, and `predicted` must be identical")

    logrank_core <- function(
        cutoff
    ){
        groups <- ifelse(predicted >= cutoff, 1, 0)
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

scores_rank <- function(predicted, actual) {
    tibble::tibble(
        rank(-predicted),
        predicted,
        ifelse(actual == 1, "high", "low")
    )
}