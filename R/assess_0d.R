#' @title Assess a model on a single data set with 0-dim
#' @description Map every model to a 0-dim metric (a single real number) and store 
#' these metrics together with some more analytics across the splits in a sorted 
#' tibble. The ground truth is time to event less or equal a cutoff.
#' @inheritParams assess_2d
#' @param ass_spec_0d AssSpec0d S3 object. The specifications for the 0-dim metric. 
#' See the constructor [`AssSpec0d()`] for details.
#' @return numeric vector. For every split index the desired metric.
#' @details The AssSpec0d S3 class is tailored for this function.
#' @export
assess_0d <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec,
    ass_spec_0d,
    quiet = FALSE,
    msg_prefix = ""
){
    prep <- prepare_and_predict(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec,
        lambda = ass_spec_0d$lambda,
        pivot_time_cutoff = ass_spec_0d$pivot_time_cutoff,
        benchmark_col = ass_spec_0d$benchmark
    )
    core <- function(i){
        predicted <- prep[["predicted"]][[i]]
        actual <- prep[["actual"]][[i]]
        pa <- intersect_by_names(predicted, actual, rm_na = TRUE)
        do.call(ass_spec_0d$metric, list(
            "predicted" =  pa[[1]],
            "actual" = pa[[2]]
            )
        )
    }
    metric <- vapply(seq_along(model_spec$split_index), core, numeric(1))
    return(metric)
}


get_auc <- function(
    predicted,
    actual
){
    pred_obj <- ROCR::prediction(predictions = predicted, labels = actual)
    res <- ROCR::performance(pred_obj, measure = "auc")
    return(res@y.values[[1]])
}