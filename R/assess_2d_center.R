#' @title Assess multiple models on a single data set in 2D plot
#' @description Assess the performance of multiple models on a single data set. 
#' For every 
#' * cohort (train and test),
#' * model and
#' * time cutoff this model specifies,
#' call [`assess_2d()`]. In addition, plot `x_metric` vs. `y_metric` attribute
#' of `AssSpec2d` for all models in a sinlge plot ("comparison plot").
#' @param model_list list of ModelSpec objects. Assess these models.
#' @param data DataSpec object. Assess on this data set.
#' @param ass2d AssSpec2d object. Specify the final comparison plot. 
#' We derive the `AssSpec2d` for the single plots in a reasonable way from it.
#' @param cohorts character vector, a subset of `c("train", "test")`. Assess on these
#' cohorts. Default is `c("train", "test")`.
#' @param model_tree_mirror character vector of lengtgh 2. This answers the question 
#' where to store the generated assessment files.
#' * For the train cohort, it is the model directory (where the fit lies as an rds 
#'  file).
#' * For the test cohort, it is the same directory as for the train cohort, with one 
#'  exception: we mirror the file path according to `model_tree_mirror`, i.e., we
#'  replace the first element of `model_tree_mirror` by the second element. Default 
#'  is `c("models", "results")`.
#' @param comparison_plot logical. Whether to generate the above mentioned comparison 
#' plot. Default is `TRUE`.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @export
assess_2d_center <- function(
    ass2d,
    model_list,
    data,
    cohorts = c("test", "train"),
    model_tree_mirror = c("models", "results"),
    comparison_plot = TRUE,
    quiet = FALSE
){
    cohorts <- match.arg(cohorts, several.ok = TRUE)

    perf_tbls <- list()

    data <- read(data)
    expr_mat <- data[["expr_mat"]]
    pheno_tbl <- data[["pheno_tbl"]]

    if(!quiet) message("\nASSESSING ON ", data$name)
    for(cohort in cohorts){
        if(!quiet) message("# On ", cohort, " cohort")
        data$cohort <- cohort
        for(model in model_list){
            if(!quiet) message("## ", model$name)
            for(time_cutoff in model$time_cutoffs){
                if(!quiet) message("### At time cutoff ", time_cutoff)
                model_cutoff <- at_time_cutoff(model, time_cutoff)
                this_as2 <- self$infer(
                    ass2d = ass2d,
                    model = model_cutoff,
                    data = data,
                    model_tree_mirror = model_tree_mirror
                )
                this_as2$assess(
                    data = data,
                    model = model_cutoff,
                    quiet = quiet,
                    msg_prefix = "#### "
                )
                perf_tbls[[model_cutoff$name]] <- this_as2$data
            }
        }
        cohort_as2 <- ass2d
        cohort_as2$data <- dplyr::bind_rows(perf_tbls)
        if(comparison_plot){
            if(cohort == "test")
                cohort_as2$file <- mirror_path(
                    filepath = ass2d$file,
                    mirror = model_tree_mirror
                )
            if(is.null(cohort_as2$title))
                cohort_as2$title <- paste0(
                    data$name, " ", data$cohort, ", ", data$time_to_event_col,
                    " < ", cohort_as2$pivot_time_cutoff
                )
            plot_2d_metric(
                ass2d = cohort_as2,
                quiet = TRUE
            )
            if(!quiet)
                message("# Saving comparative performance plot to ", cohort_as2$file)
        }
    }
}