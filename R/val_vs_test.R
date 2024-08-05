#' @title Compare validation and test error
#' @description Plot rank of validation error versus test error - validation
#' error for a list of models.
#' @param model_list list of Model objects.
#' @param data Data object.
#' @param error_fun function. Error function to compare validation and test error.
#' @param regex1 character vector. Regular-expression patterns. For every model, we assign it the 
#' group in `name1` corresponding to the first pattern its name matches and color it accordingly.
#' @param regex2 character vector. Regular-expression patterns. For every model, we assign it the
#' group in `name2` corresponding to the first pattern its name matches and shape it accordingly.
#' @param name1,name2 character vectors. The groups as explained in `regex1` and `regex2`.
#' @param legendtitle1,legendtitle2 string. Legend titles for the grouping according to `regex1` 
#' and `regex2`.
#' @param correlation_label logical. Whether to show the correlation coefficient between the
#' validation and test error as a label in the plot.
#' @param file string. File name to save the plot to.
#' @param return_type string. Either "ggplot" or "tibble". See return section for details.
#' @param colors character vector. Colors used for points.
#' @param width,height numeric. Width and height of the stored plot in inches.
#' @param plot_theme ggplot2 theme. Theme to apply to the plot. If it's `NULL`, use the ggplot2 
#' default theme.
#' @param quiet logical. Whether to suppress messages.
#' @return ggplot object if `return_type == "ggplot"` or the tibble underlying the plot if 
#' `return_type == "tibble"`.
#' @export
val_vs_test <- function(
    model_list,
    data,
    error_fun,
    regex1 = NULL,
    regex2 = NULL,
    name1 = NULL,
    name2 = NULL,
    legendtitle1 = "spot 1",
    legendtitle2 = "spot 2",
    correlation_label = TRUE,
    file = NULL,
    return_type = c("ggplot", "tibble"),
    plot_theme = ggplot2::theme_minimal(),
    colors = NULL,
    width = 7,
    height = 4,
    quiet = FALSE
) {
    return_type <- match.arg(return_type)
    n_models <- length(model_list)
    spots_regex <- list(regex1, regex2)
    spots_name <- list(name1, name2)
    for (i in seq_along(spots_regex)) {
        if (!is.null(spots_regex[[i]])) {
            spots_regex[[i]] <- c(spots_regex[[i]], ".")
            if (is.null(spots_name[[i]]))
                spots_name[[i]] <- spots_regex[[i]]
            else 
                spots_name[[i]] <- c(spots_name[[i]], "other") 
        }
    }
    tbl <- tibble::tibble(
        "validation error" = rep(NA, n_models),
        "test error" =  rep(NA, n_models),
        "spot 1" = character(n_models),
        "spot 2" = character(n_models)
    )
    names(tbl)[3:4] <- c(legendtitle1, legendtitle2)
    if (is.null(regex2)) tbl <- tbl[, -4]
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        if (!quiet) message("At model ", model$name)
        test_cohort <- data$cohort
        cohorts <- c("val_predict", test_cohort)
        skip <- FALSE
        for (j in seq_along(cohorts)) {
            data$cohort <- cohorts[j]
            yy <- model$predict(data, quiet = quiet)
            if (is.null(yy)) {
                skip <- TRUE
                break
            }
            yy <- intersect_by_names(yy[[1]], yy[[2]], rm_na = c(TRUE, TRUE))
            tbl[i, j] <- error_fun(y_hat = yy[[1]], y = yy[[2]])
        }
        if (skip) next
        idx <- 0
        for (s in seq_along(spots_regex)) {
            if (is.null(spots_regex[[s]])) next
            for (j in seq_along(spots_regex[[s]])) {
                if (stringr::str_detect(model$name, spots_regex[[s]][j])) {
                    tbl[i, 2+s] <- spots_name[[s]][j]
                    break
                }
            }
        }
    }
    tbl <- tbl[!is.na(tbl[[1]]), ] # Remove skipped rows
    if (return_type == "tibble") return(tbl)

    showtext::showtext_auto()
    plt <- ggplot2::ggplot(data = tbl, mapping = ggplot2::aes(
        x = .data[["validation error"]],
        y = .data[["test error"]]
    )) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgray") +
        plot_theme +
        ggplot2::guides(
            color = ggplot2::guide_legend(byrow = TRUE),
            shape = ggplot2::guide_legend(byrow = TRUE)
        ) +
        ggplot2::scale_shape_manual(values = c(4, 2, 5, 6, 0, 3))
    if (!is.null(regex1)) plt <- plt + ggplot2::aes(color = .data[[legendtitle1]])
    if (!is.null(regex2)) 
        plt <- plt + ggplot2::aes(shape = .data[[legendtitle2]]) + 
            ggplot2::geom_point()
    else 
        plt <- plt + ggplot2::geom_point(shape = 4)
    if (correlation_label) 
        plt <- plt +
        ggplot2::geom_label(
            x = max(tbl[["validation error"]]),
            y = max(tbl[["test error"]]),
            label = paste0("rho = ", round(stats::cor(tbl[["validation error"]], 
                tbl[["test error"]]), 2)),
            inherit.aes = FALSE,
            hjust = 1,
            vjust = 1,
            alpha = 1,
            fill = "white",
            family = plot_theme$text$family
        )
    if (!is.null(colors)) plt <- plt + ggplot2::scale_color_manual(values = colors)
    if (!is.null(file)) ggplot2::ggsave(file, plt, width = width, height = height, dpi = 300)
    showtext::showtext_auto(FALSE)
    invisible(plt)
}