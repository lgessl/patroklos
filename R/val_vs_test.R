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
    plot_theme = NULL,
    quiet = FALSE
) {
    if (is.null(spotlight_name)) spotlight_name <- spotlight_regex
    n_models <- length(model_list)
    tbl <- tibble::tibble(
        val_error = numeric(n_models),
        test_error = numeric(n_models),
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
    tbl[["val error - test error"]] <- tbl[["val_error"]] - tbl[["test_error"]]
    tbl[["val error rank"]] <- rank(tbl[["val_error"]])

    plt <- ggplot2::ggplot(data = tbl, mapping = ggplot2::aes(
        x = .data[["val error rank"]],
        y = .data[["val error - test error"]],
        color = .data[["spotlight"]]
    )) + ggplot2::geom_point()
    if (!is.null(plot_theme)) plt <- plt + plot_theme
    ggplot2::ggsave(file, plt, width = 7, height = 4, dpi = 300)
    invisible(plt)
}