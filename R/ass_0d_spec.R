Ass0dSpec <- function(
    metric,
    pivot_time_cutoff,
    lambda = "lambda.min",
    benchmark = NULL
){
    stopifnot(is.character(metric))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(pivot_time_cutoff > 0)
    stopifnot(is.character(lambda) || is.numeric(lambda))
    stopifnot(is.null(benchmark) || is.character(benchmark))

    ass_0d_spec <- list(
        "metric" = metric,
        "pivot_time_cutoff" = pivot_time_cutoff,
        "lambda" = lambda,
        "benchmark" = benchmark
    )
    return(structure(ass_0d_spec, class = "Ass0dSpec"))
}