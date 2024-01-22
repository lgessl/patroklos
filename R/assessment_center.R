#' @title Assess multiple models on a single data set
#' @description Assess the performance of multiple models on a single data set. 
#' For every 
#' * cohort (train and test),
#' * model and
#' * time cutoff this model specifies,
#' call [`assess_model()`]. In addition, plot `x_metric` vs. `y_metric` attribute
#' of `PerfPlotSpec` for all models in a sinlge plot ("comparison plot").
#' @param model_spec_list list of ModelSpec objects. Assess these models.
#' @param data_spec DataSpec object. Assess on this data set.
#' @param perf_plot_spec PerfPlotSpec object. Specify the final comparison plot. 
#' We derive the `PerfPlotSpec` for the single plots in a reasonable way from it.
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
assessment_center <- function(
    model_spec_list,
    data_spec,
    perf_plot_spec,
    cohorts = c("train", "test"),
    model_tree_mirror = c("models", "results"),
    comparison_plot = TRUE,
    quiet = FALSE
){
    cohorts <- match.arg(cohorts, several.ok = TRUE)

    perf_tbls <- list()

    data <- read(data_spec)
    expr_mat <- data[["expr_mat"]]
    pheno_tbl <- data[["pheno_tbl"]]
    perf_plot_spec$model_tree_mirror <- model_tree_mirror

    message("\nASSESSING ON ", data_spec$name)
    for(cohort in c("train", "test")){
        if(!quiet) message("# On ", cohort, " cohort")
        data_spec$cohort <- cohort
        cohort_pps <- perf_plot_spec
        for(model_spec in model_spec_list){
            if(!quiet) message("## ", model_spec$name)
            for(time_cutoff in model_spec$time_cutoffs){
                if(!quiet) message("### At time cutoff ", time_cutoff)
                ms_cutoff <- at_time_cutoff(model_spec, time_cutoff)
                this_pps <- infer_pps(
                    perf_plot_spec = cohort_pps,
                    model_spec = ms_cutoff,
                    data_spec = data_spec
                )
                this_pps <- assess_model(
                    expr_mat = expr_mat,
                    pheno_tbl = pheno_tbl,
                    data_spec = data_spec,
                    model_spec = ms_cutoff,
                    perf_plot_spec = this_pps,
                    quiet = quiet,
                    msg_prefix = "#### "
                )
                perf_tbls[[ms_cutoff$name]] <- this_pps$data
            }
        }
        cohort_pps$data <- dplyr::bind_rows(perf_tbls)
        if(cohort == "test")
            cohort_pps$fname <- mirror_directory(
                filepath = cohort_pps$fname,
                mirror = cohort_pps$model_tree_mirror
            )
        if(is.null(cohort_pps$title))
            cohort_pps$title <- paste0(
                data_spec$name, " ", data_spec$cohort, ", ", data_spec$time_to_event_col,
                " < ", cohort_pps$pivot_time_cutoff
            )
        if(comparison_plot){
            plot_2d_metric(
                perf_plot_spec = cohort_pps,
                quiet = TRUE
            )
            if(!quiet)
                message("# Saving comparative performance plot to ", cohort_pps$fname)
        }

    }

}