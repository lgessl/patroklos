#' @title An S3 class specifying how to assess with a 0-dim metric (like AUC)
#' @inheritParams AssSpec2d
#' @param metric character. Name of the function used to calculate the metric. It must 
#' take two numeric vectors of the same length as arguments: 
#' 1. `predicted`, the predictions by a model, higher values indicate a more positive 
#' prediction,
#' 2. `actual`, the true values with the positive class coded as 1, the negative one as 0. 
#' @param file string. The name of the csv file to save the results to. If `NULL`,
#' the results are not saved. If you supply an AssSpec0d object to `assess_0d_center()`,
#' specify the path for the *train* cohort.
#' @param round_digits numeric. The number of digits to round the results to. Default is `3`.
#' @return An AssSpec0d S3 object.
#' @export
AssSpec0d <- function(
    metric,
    pivot_time_cutoff,
    lambda = "lambda.min",
    benchmark = NULL,
    file = NULL,
    round_digits = 3
){
    stopifnot(is.character(metric))
    stopifnot(is.numeric(pivot_time_cutoff))
    stopifnot(pivot_time_cutoff > 0)
    stopifnot(is.character(lambda) || is.numeric(lambda))
    stopifnot(is.null(benchmark) || is.character(benchmark))
    stopifnot(is.null(file) || is.character(file))
    stopifnot(is.numeric(round_digits) && round_digits >= 0)

    ass_spec_0d <- list(
        "metric" = metric,
        "pivot_time_cutoff" = pivot_time_cutoff,
        "lambda" = lambda,
        "benchmark" = benchmark,
        "file" = file,
        "round_digits" = round_digits
    )
    return(structure(ass_spec_0d, class = "AssSpec0d"))
}