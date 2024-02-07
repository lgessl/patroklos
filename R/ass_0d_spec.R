AssSpec0d <- function(
    metric,
    pivot_time_cutoff,
    lambda = "lambda.min",
    benchmark = NULL,
    file = NULL
){
    stopifnot(is.character(metric))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(pivot_time_cutoff > 0)
    stopifnot(is.character(lambda) || is.numeric(lambda))
    stopifnot(is.null(benchmark) || is.character(benchmark))
    stopifnot(is.null(file) || is.character(file))

    ass_spec_0d <- list(
        "metric" = metric,
        "pivot_time_cutoff" = pivot_time_cutoff,
        "lambda" = lambda,
        "benchmark" = benchmark,
        "file" = file
    )
    return(structure(ass_spec_0d, class = "AssSpec0d"))
}