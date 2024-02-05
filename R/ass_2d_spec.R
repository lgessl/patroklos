new_Ass2dSpec <- function(
    fname,
    x_metric,
    y_metric,
    pivot_time_cutoff,
    lambda,
    benchmark,
    ci_level,
    fellow_csv,
    scores_plot,
    show_plots,
    title,
    x_lab,
    y_lab,
    xlim,
    ylim,
    smooth_method,
    smooth_benchmark,
    scale_x,
    scale_y,
    hline,
    vline,
    text,
    alpha,
    colors,
    width,
    height,
    units
){
    stopifnot(is.character(fname))
    stopifnot(is.character(x_metric))
    stopifnot(is.character(y_metric))
    stopifnot(is.null(pivot_time_cutoff) || is.numeric(pivot_time_cutoff))
    stopifnot(is.character(lambda) || is.numeric(lambda))
    stopifnot(is.character(benchmark) || is.null(benchmark))
    stopifnot(is.numeric(ci_level) && ci_level >= 0 && ci_level <= 1)
    stopifnot(is.logical(fellow_csv))
    stopifnot(is.logical(scores_plot))
    stopifnot(is.logical(show_plots))
    stopifnot(is.character(title) || is.null(title))
    stopifnot(is.character(x_lab))
    stopifnot(is.character(y_lab))
    stopifnot(is.numeric(xlim) || is.null(xlim))
    stopifnot(is.numeric(ylim) || is.null(ylim))
    stopifnot(is.character(smooth_method) || is.null(smooth_method) ||
        is.function(smooth_method))
    stopifnot(is.logical(smooth_benchmark))
    stopifnot(is.null(hline) || is.list(hline))
    stopifnot(is.null(vline) || is.list(vline))
    stopifnot(is.null(text) || is.list(text))
    stopifnot(is.numeric(alpha) && alpha >= 0 && alpha <= 1)
    stopifnot(is.character(colors) || is.null(colors))
    stopifnot(is.numeric(width) && width > 0)
    stopifnot(is.numeric(height) && height > 0)
    stopifnot(is.character(units))

    ass_2d_spec_list <- list(
        "fname" = fname,
        "x_metric" = x_metric,
        "y_metric" = y_metric,
        "pivot_time_cutoff" = pivot_time_cutoff,
        "lambda" = lambda,
        "benchmark" = benchmark,
        "ci_level" = ci_level,
        "fellow_csv" = fellow_csv,
        "scores_plot" = scores_plot,
        "show_plots" = show_plots,
        "title" = title,
        "x_lab" = x_lab,
        "y_lab" = y_lab,
        "xlim" = xlim,
        "ylim" = ylim,
        "smooth_method" = smooth_method,
        "smooth_benchmark" = smooth_benchmark,
        "scale_x" = scale_x,
        "scale_y" = scale_y,
        "hline" = hline,
        "vline" = vline,
        "text" = text,
        "alpha" = alpha,
        "colors" = colors,
        "width" = width,
        "height" = height,
        "units" = units
    )
    return(structure(ass_2d_spec_list, class = "Ass2dSpec"))
}

#' @title Create a Ass2dSpec object
#' @description A Ass2dSpec object holds all the info on how to assess how well a 
#' model filters high-risk patients. Its base object is a list. The core of the 
#' assessment is a scatter plot of two performance measures. Moreover, a Ass2dSpec
#' object holds the commands (usually bools) on whether do more assessments (that 
#' require considerably less specifications).
#' @param fname string. The name of the file to save the plot to.
#' @param x_metric string. The name of the performance measure to be plotted on the x-axis.
#' All measures that can be passed to the `x.measure` parameter of [ROCR::performance()] are
#' valid.
#' @param y_metric string. The name of the performance measure to be plotted on the y-axis.
#' * All measures that can be passed to the `measure` parameter of [ROCR::performance()] are
#' valid.
#' * `"logrank"` for the p-values of a logrank test. In this case, `x_metric` must be
#' `"prevalence"` or `"rpp"`.
#' @param pivot_time_cutoff numeric. Time-to-event threshold that divides samples into a
#' high/low-risk (time to event below/above `pivot_time_cutoff`) group. Model performance
#' will be measured in terms of how well it can separate these two groups as a binary
#' classifier. Default is NULL, in which case the `pivot_time_cutoff` of an accompanying
#' DataSpec is used.
#' @param lambda numeric or string. The lambda regularization parameter of the model to
#' evaluate. See [prepare_and_predict()] for details. Default is `"lambda.min"`.
#' @param benchmark string. Column in the test pheno holding the numeric benchmark values.
#' Default is `"ipi"` (international prognostic index for DLBCL).
#' @param ci_level numeric in \[0, 1\]. The level used to calculate confidence intervals.
#' Default is `0.95`.
#' @param fellow_csv logical. If passed to [assess_2d_center()], whether to also create a
#' csv file for every model-data pair. Default is `TRUE`.
#' @param scores_plot logical. Display the ordered scores output by the model in a scatter 
#' plot. Default is `TRUE`.  
#' @param show_plots logical. Whether to show the plots after creation in an interactive 
#' graphics device. Default is `FALSE`. 
#' @param title string. The title of the plot. Default is `NULL`.
#' @param x_lab,y_lab string. Axis labels. Default is `x_metric` and `y_metric`, respectively.
#' @param xlim,ylim numeric vector of length 2. The limits for both axes. Default is
#' `c(-Inf, Inf)`, i.e. no contraints.
#' @param smooth_method string or function. Smooth method to plot an additional smoothed graph.
#' If `NULL`, no smoothing. Else we pass `smooth_method` as the `method` parameter to
#' [ggplot2::geom_smooth()]. Default is `"loess"`.
#' @param smooth_benchmark logical. Whether to also smooth the benchmark data. Default is
#' `FALSE`.
#' @param scale_x,scale_y string or transformation object (see [`scales::trans_new`] for the 
#' latter). The scale of the axes, we will pass them to the `trans` paramter of 
#' [`ggplot2::scale_x_continuous()`], [`ggplot2::scale_y_continuous()], respectively. 
#' Default is `"identity"`.
#' @param vline,hline list or NULL. Vertical/horizontal lines to be added to the plot. A list
#' holding the arguments to pass to [`ggplot2::geom_vline()`] and [`ggplot2::geom_hline()`],
#' respectively. Default is `NULL`.
#' @param text list or NULL. Text label added to the plot. A list holding the arguments to
#' pass to [`ggplot2::geom_label()`]. Default is `NULL`. 
#' @param alpha numeric in \[0, 1\]. The alpha value for the points and lines in the 
#' plot.
#' @param width numeric. The width of the plot in `units`. Default is `7`.
#' @param height numeric. The height of the plot in `units`. Default is `4`.
#' @param colors character vector. The colors to be used for the different models.
#' Default is `NULL`, which means that the default colors of `ggplot2` will be used.
#' @param units string. The units of `width` and `height`. Default is `"in"` (inches).
#' @return A Ass2dSpec S3 object.
#' @export
Ass2dSpec <- function(
    fname,
    x_metric,
    y_metric,
    pivot_time_cutoff = NULL,
    lambda = "lambda.min",
    benchmark = NULL,
    ci_level = 0.95,
    fellow_csv = FALSE,
    scores_plot = FALSE,
    show_plots = FALSE,
    title = NULL,
    x_lab = NULL,
    y_lab = NULL,
    xlim = c(-Inf, Inf),
    ylim = c(-Inf, Inf),
    smooth_method = NULL,
    smooth_benchmark = FALSE,
    scale_x = "identity",
    scale_y = "identity",
    vline = NULL,
    hline = NULL,
    text = NULL,
    alpha = 0.5,
    width = 7,
    height = 4,
    colors = NULL,
    units = "in"
){
    if(is.null(x_lab)){
        x_lab <- x_metric
    }
    if(is.null(y_lab)){
        y_lab <- y_metric
    }
    if(y_metric == "logrank"){
        if(!any(stringr::str_detect(c("prevalence", "rpp"), x_metric)))
            stop("For `y_metric` = 'logrank', `x_metric` must be 'prevalence' or 'rpp'.")
    }
    ass_2d_spec <- new_Ass2dSpec(
        fname = fname,
        x_metric = x_metric,
        y_metric = y_metric,
        pivot_time_cutoff = pivot_time_cutoff,
        lambda = lambda,
        benchmark = benchmark,
        ci_level = ci_level,
        fellow_csv = fellow_csv,
        scores_plot = scores_plot,
        x_lab = x_lab,
        y_lab = y_lab,
        xlim = xlim,
        ylim = ylim,
        smooth_method = smooth_method,
        smooth_benchmark = smooth_benchmark,
        scale_x = scale_x,
        scale_y = scale_y,
        vline = vline,
        hline = hline,
        text = text,
        alpha = alpha,
        width = width,
        height = height,
        units = units,
        title = title,
        colors = colors,
        show_plots = show_plots
    )
    return(ass_2d_spec)
}


infer_pps <- function(
    ass_2d_spec,
    model_spec,
    data_spec
){
    # Prepare for assess_2d()
    this_pps <- ass_2d_spec
    this_pps$fname <- file.path(
        model_spec$directory,
        paste0(
            ass_2d_spec$x_metric, "_vs_", ass_2d_spec$y_metric,
            stringr::str_extract(ass_2d_spec$fname, "\\..+$")
        )
    )
    if(data_spec$cohort == "test")
        this_pps$fname <- mirror_directory(
            filepath = this_pps$fname,
            mirror = ass_2d_spec$model_tree_mirror
        )
    this_pps$title <- paste0(
        data_spec$name, " ", data_spec$cohort, ", ",
        data_spec$time_to_event_col, " < ", this_pps$pivot_time_cutoff
    )
    return(this_pps)
}