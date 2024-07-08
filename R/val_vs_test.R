#' @title Compare validation and test error
#' @description Plot rank of validation error versus test error - validation
#' error for a list of models.
#' @param model_list list of Model objects.
#' @param data Data object.
#' @param error_fun function. Error function to compare validation and test error.
#' @param spotlight_regex string. Regular expression to identify models to
#' highlight.
#' @param file string. File name to save the plot to.
#' @param spotlight_name string. Name to give to the models identified by
#' `spotlight_regex` in legend.
#' @param colors character vector. Colors used for points.
#' @param width,height numeric. Width and height of the plot in inches.
#' @param plot_theme ggplot2 theme. Theme to apply to the plot.
#' @param quiet logical. Whether to suppress messages.
#' @return ggplot2 object.
#' @export
val_vs_test <- function(
    model_list,
    data,
    error_fun,
    spotlight_regex,
    file,
    spotlight_name = NULL,
    plot_theme = ggplot2::theme_minimal(),
    colors = NULL,
    width = 7,
    height = 4,
    quiet = FALSE
) {
    if (is.null(spotlight_name)) spotlight_name <- spotlight_regex
    n_models <- length(model_list)
    tbl <- tibble::tibble(
        "validation error" = numeric(n_models),
        "test error" = numeric(n_models),
        spotlight = character(n_models)
    )
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        model <- readRDS(file.path(model$directory, model$file))
        test_cohort <- data$cohort
        cohorts <- c("val_predict", test_cohort)
        for (j in seq_along(cohorts)) {
            data$cohort <- cohorts[j]
            yy <- model$predict(data, quiet = quiet)
            yy <- intersect_by_names(yy[[1]], yy[[2]], rm_na = c(TRUE, TRUE))
            tbl[i, j] <- error_fun(y_hat = yy[[1]], y = yy[[2]])
        }
        tbl[i, "spotlight"] <- ifelse(
            stringr::str_detect(model$name, spotlight_regex),
            spotlight_name,
            "other"
        )
    }

    if (length(unique(tbl[["spotlight"]])) > 1)
        plt <- ggplot2::ggplot(data = tbl, mapping = ggplot2::aes(
            x = .data[["validation error"]],
            y = .data[["test error"]],
            color = .data[["spotlight"]]
        )) + ggplot2::geom_point()
    else {
        if (is.null(colors)) colors <- "black"
        plt <- ggplot2::ggplot(data = tbl, mapping = ggplot2::aes(
            x = .data[["validation error"]],
            y = .data[["test error"]]
        )) + ggplot2::geom_point(color = colors[1])
    }
    plt <- plt +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    ggplot2::geom_label(
        x = max(tbl[["validation error"]]),
        y = max(tbl[["test error"]]),
        label = paste("rho =", round(stats::cor(tbl[["validation error"]], 
            tbl[["test error"]]), 2)),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = 1,
        alpha = 1,
        fill = "white",
        family = plot_theme$text$family
    ) + plot_theme
    if (!is.null(colors)) plt <- plt + ggplot2::scale_color_manual(values = colors)
    ggplot2::ggsave(file, plt, width = width, height = height, dpi = 300)
    invisible(plt)
}