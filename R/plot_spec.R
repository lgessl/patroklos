new_PlotSpec <- function(
    fname,
    x_axis,
    y_axis,
    width,
    height,
    units,
    title,
    colors,
    show_plots
){
    stopifnot(is.character(fname))
    stopifnot(is.character(x_axis))
    stopifnot(is.character(y_axis))
    stopifnot(is.numeric(width))
    stopifnot(is.numeric(height))
    stopifnot(is.character(units))
    stopifnot(is.character(title) || is.null(title))
    stopifnot(is.character(colors) || is.null(colors))
    stopifnot(is.logical(show_plots))

    plot_spec_list <- list(
        "fname" = fname,
        "x_axis" = x_axis,
        "y_axis" = y_axis,
        "width" = width,
        "height" = height,
        "units" = units,
        "title" = title,
        "colors" = colors
    )
    return(structure(plot_spec_list, class = "PlotSpec"))
}

PlotSpec <- function(
    fname,
    x_axis = c("prevalence", "true positive rate"),
    y_axis = c("precision", "false positive rate"),
    width = 7,
    height = 4,
    units = "in",
    title = NULL,
    colors = NULL,
    show_plots = FALSE
){
    x_axis <- match.arg(x_axis)
    y_axis <- match.arg(y_axis)
    plot_spec <- new_PlotSpec(
        fname,
        x_axis,
        y_axis,
        width,
        height,
        units,
        title,
        colors,
        show_plots
    )
    return(plot_spec)
}