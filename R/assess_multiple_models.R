#' @title Assess multiple models and (optionally) compare them in a single 
#' plot
#' @description Assess the performance of multiple models and (optionally)
#' multiple data sets, finally and optionally compare them in a single plot. 
#' Consider the model a classifier for high-risk versus low-risk patients 
#' (depending on `per_plot_spec$pfs_leq`).
#' For every model-data pair, plot one classifier property on the x-axis 
#' against another one on the y-axis. If all models are evaluated on the same
#' data set, the comparative plot can also contain a benchmark classifier given 
#' via a score in the pheno data (e.g. the International Prognostic Index (IPI) 
#' for DLBCL).
#' @param model_spec_list list of ModelSpec objects. The models to be assessed.
#' @param data_spec_list list of DataSpec objects. Data sets corresponing to the
#' models in `model_spec_list`. If is has length 1, the same data will be used 
#' to evaluate every model specified in `model_spec_list`. If more than one data 
#' set is given, the benchmark classifier will not be plotted to avoid an 
#' overcrowded plot.
#' @param perf_plot_spec PerfPlotSpec object. The specifications for the single
#' and the comparison plot.
#' @param single_plots logical. Whether to also generate individual plots for
#' every model-data pair. Default is `TRUE`.
#' @param model_tree_mirror character vector of lengtgh 2. For the individual plots,
#' mirror the model file tree by replacing the directory `tree_mirror[1]` in the
#' model file path by `tree_mirror[2]`. Default is `c("models", "results")`.
#' @param comparison_plot logical. Whether to compare all models in a single plot.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @return A tibble containing the performance measures for all models and the
#' benchmark classifier (if given and included in the plot).
#' @details If `single_plots == TRUE`, create a single plot for every model-data 
#' pair. For every such plot, a PlotSpec object inheriting
#' from `perf_plot_spec` is derived in a reasonable way. If 
#' `perf_plot_spec$fellow_csv` is `TRUE`, a csv table holding the plot the data 
#' will show up in the same directory under the same name except
#' for the file extension (`.csv` instead of `.pdf`). Use `perf_plot_spec$fname` and
#' `perf_plot_spec$title` for the comparison plot.
#' @export
assess_multiple_models <- function(
    model_spec_list,
    data_spec_list,
    perf_plot_spec,
    single_plots = TRUE,
    model_tree_mirror = c("models", "results"),
    comparison_plot = TRUE,
    quiet = FALSE
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
    data_spec <- NULL
    expr_mat <- NULL
    pheno_tbl <- NULL
    for(i in seq_along(model_spec_list)){
        # Extract and read in
        if(!one_data || (one_data && i == 1L)){
            data_spec <- data_spec_list[[i]]
            data <- read(data_spec)
            expr_mat <- data[["expr_mat"]]
            pheno_tbl <- data[["pheno_tbl"]]
        }
        model_spec <- model_spec_list[[i]]

        # Prepare for assess_model()
        single_pps <- perf_plot_spec
        if(single_plots){
            dir <- mirror_directory(
                filepath = model_spec$save_dir,
                mirror = model_tree_mirror
            )
            if(!dir.exists(dir))
                dir.create(dir, recursive = TRUE)
            base_fname <- stringr::str_c(
                perf_plot_spec$x_metric, "_vs_", perf_plot_spec$y_metric, ".pdf"
            )
            single_pps$fname <- file.path(dir, base_fname)
            single_pps$title <- stringr::str_c(
                data_spec$name, ", pfs < ", model_spec$pfs_leq
            )
        }

        single_tbls <- assess_model(
            expr_mat = expr_mat,
            pheno_tbl = pheno_tbl,
            data_spec = data_spec,
            model_spec = model_spec,
            perf_plot_spec = single_pps,
            plots = single_plots,
            quiet = quiet
        )
        perf_tbls[names(single_tbls)] <- single_tbls

        if(!one_data) # Remove if more than one data set is given
            perf_tbls[[perf_plot_spec$benchmark]] <- NULL
    }

    perf_tbl <- dplyr::bind_rows(perf_tbls, .id = "model")
    perf_plot_spec$data <- perf_tbl
    if(comparison_plot){
        plot_perf_metric(
            perf_plot_spec = perf_plot_spec,
            quiet = TRUE
        )
        message("Saving comparative performance plot to ", perf_plot_spec$fname)
    }

    invisible(perf_tbl)
}


#' @title Call [assess_multiple_models()] on training and test data
#' @description Given data specifications for training and test data, a list
#' of model specifications and a plot specification for the training data, call 
#' `assess_multiple_models()` twice: for the training data and, after mirroring
#' the comparison plot file path to the results directory, for the test data.
#' @inheritParams assess_multiple_models
#' @param data_spec_train DataSpec object for the training data.
#' @param data_spec_test DataSpec object for the test data.
#' @param perf_plot_spec_train PerfPlotSpec object for the training data. The same
#' PerfPlotSpec object will be used for the test data, except that the file path
#' will be mirrored to the results directory.
#' @param model_tree_mirror character vector of lengtgh 2. See [assess_multiple_models()].
#' Here, we also use it to mirror the `fname` attribute of `perf_plot_spec_train` to
#' obtain the `PerfPlotSpec` object for the test data.
#' @export
assess_train_and_test <- function(
    model_spec_list,
    data_spec_train,
    data_spec_test,
    perf_plot_spec_train,
    model_tree_mirror = c("models", "results"),
    single_plots = TRUE,
    comparison_plot = TRUE,
    quiet = FALSE
){
    # Train
    if(!quiet)
        message("Calling assess_multiple_models() on training data")
    assess_multiple_models(
        model_spec_list = model_spec_list,
        data_spec_list = list(data_spec_train),
        perf_plot_spec = perf_plot_spec_train,
        single_plots = single_plots,
        model_tree_mirror = rep(model_tree_mirror[1], 2),
        comparison_plot = comparison_plot,
        quiet = quiet
    )

    # Test
    perf_plot_spec_test <- perf_plot_spec_train
    perf_plot_spec_test$fname <- mirror_directory(
        filepath = perf_plot_spec_test$fname,
        mirror = model_tree_mirror
    )
    if(!quiet)
        message("Calling assess_multiple_models() on test data")
    assess_multiple_models(
        model_spec_list = model_spec_list,
        data_spec_list = list(data_spec_test),
        perf_plot_spec = perf_plot_spec_test,
        single_plots = single_plots,
        model_tree_mirror = model_tree_mirror,
        comparison_plot = comparison_plot,
        quiet = quiet
    )
}