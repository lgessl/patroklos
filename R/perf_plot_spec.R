new_PerfPlotSpec <- function(
    fname,
    x_metric,
    y_metric,
    benchmark,
    model_names_in_legend,
    also_single_plots,
    single_csvs,
    x_lab,
    y_lab,
    width,
    height,
    units,
    title,
    colors,
    show_plots
){
    stopifnot(is.character(fname))
    stopifnot(is.character(x_metric))
    stopifnot(is.character(y_metric))
    stopifnot(is.character(benchmark))
    stopifnot(is.logical(also_single_plots))
    stopifnot(is.numeric(width))
    stopifnot(is.numeric(height))
    stopifnot(is.character(units))
    stopifnot(is.character(title) || is.null(title))
    stopifnot(is.character(colors) || is.null(colors))
    stopifnot(is.logical(show_plots))

    perf_plot_spec_list <- list(
        "fname" = fname,
        "x_metric" = x_metric,
        "y_metric" = y_metric,
        "benchmark" = benchmark,
        "model_names_in_legend" = model_names_in_legend,
        "also_single_plots" = also_single_plots,
        "single_csvs" = single_csvs,
        "x_lab" = x_lab,
        "y_lab" = y_lab,
        "width" = width,
        "height" = height,
        "units" = units,
        "title" = title,
        "colors" = colors,
        "show_plots" = show_plots
    )
    return(structure(perf_plot_spec_list, class = "PerfPlotSpec"))
}

#' @title Create a PlotSpec object
#' @description A PlotSpec object holds all the info on how to create a model performance
#' plot and what to do with it. Its base object is a list containing:
#' @param fname string. The name of the file to save the plot to.
#' @param x_metric string. The name of the performance measure to be plotted on the x-axis.
#' All measures that can be passed to the `x.measure` parameter of `ROCR::performance()` are
#' valid.
#' @param y_metric string. The name of the performance measure to be plotted on the y-axis.
#' All measures that can be passed to the `measure` parameter of `ROCR::performance()` are
#' valid.
#' @param benchmark string. Column in the test pheno holding the numeric benchmark values.
#' Default is `"ipi"` (international prognostic index for DLBCL).
#' @param also_single_plots logical. If passed to `compare_models()`, whether to also create
#' a single plot for every model-data pair. Default is `TRUE`.
#' @param width numeric. The width of the plot in `units`. Default is `7`.
#' @param height numeric. The height of the plot in `units`. Default is `4`.
#' @param units string. The units of `width` and `height`. Default is `"in"` (inches).
#' @param title string. The title of the plot. Default is `NULL`.
#' @param colors character vector. The colors to be used for the different models.
#' Default is `NULL`, which means that the default colors of `ggplot2` will be used.
#' @param show_plots logical. Whether to show the plots after creation in an interactive 
#' graphics device. Default is `FALSE`. 
#' @return A PlotSpec S3 object.
PerfPlotSpec <- function(
    fname,
    x_metric,
    y_metric,
    benchmark = "ipi",
    model_names_in_legend = NULL,
    also_single_plots = TRUE,
    single_csvs = TRUE,
    x_lab = NULL,
    y_lab = NULL,
    width = 7,
    height = 4,
    units = "in",
    title = NULL,
    colors = NULL,
    show_plots = FALSE
){
    if(is.null(x_lab)){
        x_lab <- x_metric
    }
    if(is.null(y_lab)){
        y_lab <- y_metric
    }
    per_plot_spec <- new_PerfPlotSpec(
        fname = fname,
        x_metric = x_metric,
        y_metric = y_metric,
        benchmark = benchmark,
        model_names_in_legend = model_names_in_legend,
        also_single_plots = also_single_plots,
        single_csvs = single_csvs,
        x_lab = x_lab,
        y_lab = y_lab,
        width = width,
        height = height,
        units = units,
        title = title,
        colors = colors,
        show_plots = show_plots
    )
    return(per_plot_spec)
}