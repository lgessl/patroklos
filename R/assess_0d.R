assess_0d <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec,
    ass_0d_spec,
    quiet = FALSE,
    msg_prefix = ""
){
    prep <- prepare_and_predict(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec,
        lambda = ass_0d_spec$lambda,
        pivot_time_cutoff = ass_0d_spec$pivot_time_cutoff,
        benchmark_col = ass_0d_spec$benchmark
    )
    core <- function(i){
        predicted <- prep[["predicted"]][[i]]
        actual <- prep[["actual"]][[i]]
        available <- !is.na(predicted) & !is.na(actual)
        do.call(ass_0d_spec$metric, list(
            "predicted" =  predicted[available],
            "actual" = actual[available]
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