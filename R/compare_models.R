#' @title Compare multiple models in a single plot
#' @description Compare the performance of multiple models and (optionally)
#' multiple data sets in a single plot. Consider the model a classifier for
#' high-risk versus low-risk patients (depening on `per_plot_spec$pfs_leq`).
#' For every model-data pair, plot one classifier property on the x-axis 
#' against another one on the y-axis. If all models are evaluated on the same
#' data set, the plot can also contain a benchmark classifier given via a score
#' in the pheno data (e.g. the International Prognostic Index (IPI) for DLBCL).
#' @param model_spec_list list of ModelSpec objects. The models to be compared.
#' @param data_spec_list list of DataSpec objects. Data sets corresponing to the
#' models in `model_spec_list`. If is has length 1, the same data will be used 
#' to evaluate every model specified in `model_spec_list`. If more than one data 
#' set is given, the benchmark classifier will not be plotted to avoid an 
#' overcrowded plot.
#' @param perf_plot_spec PerfPlotSpec object. The specifications for the plot.
#' @return A tibble containing the performance measures for all models and the
#' benchmark classifier (if given and included in the plot).
#' @details If `perf_plot_spec$also_single_plots` is `TRUE`, create a single plot
#' for every model-data pair. For every such plot, a PlotSpec object inheriting
#' from `perf_plot_spec` is derived in a reasonable way and stored in the model 
#' directory under `paste0(perf_plot_spec$x_metric, _vs_, perf_plot_spec$y_metric, 
#' .pdf)`. If `perf_plot_spec$single_csvs` is `TRUE`, a csv table holding the 
#' plot the data will show up in the same directory under the same name except
#' for the file extension (`.csv` instead of `.pdf`).
compare_models <- function(
    model_spec_list,
    data_spec_list,
    perf_plot_spec
){
    perf_tbls <- list()
    # If more than one data_spec (i.e., several data sets) are given, do not plot
    # the benchmark (too crowdy) 
    one_data <- FALSE
    if(length(data_spec_list) == 1L){
        data_spec_list <- rep(data_spec_list, length(model_spec_list))
        one_data <- TRUE
    }

    # Collect tibbles, (optionally) generate individual plots, csv files
    for(i in seq_along(model_spec_list)){

        # Extract and read in
        data_spec <- data_spec_list[[i]]
        model_spec <- model_spec_list[[i]]
        data <- read(data_spec)
        expr_mat <- data[["expr_mat"]]
        pheno_tbl <- data[["pheno_tbl"]]

        # Prepare for assess_model()
        single_perf_plot_spec <- perf_plot_spec
        single_perf_plot_spec$fname <- file.path(
            model_spec$save_dir,
            stringr::str_c(perf_plot_spec$x_metric, "_vs_", 
                perf_plot_spec$y_metric, ".pdf")
        )

        single_tbls <- assess_model(
            expr_mat = expr_mat,
            pheno_tbl = pheno_tbl,
            data_spec = data_spec,
            model_spec = model_spec,
            perf_plot_spec = single_perf_plot_spec
        )
        perf_tbls[names(single_tbls)] <- single_tbls

        if(!one_data) # Remove if more than one data set is given
            perf_tbls[[perf_plot_spec$benchmark]] <- NULL
    }

    perf_tbl <- dplyr::bind_rows(perf_tbls, .id = "model")
    plot_perf_metric(
        perf_tbl = perf_tbl,
        perf_plot_spec = perf_plot_spec
    )

    return(perf_tbl)
}