#' @importFrom rlang .data
plot_2d_metric <- function(
    ass2d,
    quiet = FALSE,
    msg_prefix = ""
){
    # Extract
    data <- ass2d$data |> dplyr::distinct()
    # For legend order: bechmark last
    if(!is.null(ass2d$benchmark)){
        model_levels <- as.character(unique(data[["model"]]))
        model_levels <- sort(model_levels[model_levels != ass2d$benchmark])
        model_levels <- c(ass2d$benchmark, model_levels)
        data[["model"]] <- factor(data[["model"]], levels = model_levels)
    }
    x_metric <- ass2d$x_metric
    y_metric <- ass2d$y_metric
    # Constraint data to range
    data <- data[
        data[[x_metric]] >= ass2d$xlim[1] &
        data[[x_metric]] <= ass2d$xlim[2],
    ]
    data <- data[
        data[[y_metric]] >= ass2d$ylim[1] &
        data[[y_metric]] <= ass2d$ylim[2],
    ]
    all_data <- data
    bm_data <- NULL
    if(!is.null(ass2d$benchmark)){
        bm_data <- data[data[["model"]] == ass2d$benchmark, ]
        data <- data[data[["model"]] != ass2d$benchmark, ]
    }
    # One font family for all text
    font_family <- ass2d$theme$text$family
    if(is.null(font_family)) font_family <- ""
 
    plt <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = .data[[x_metric]], 
            y = .data[[y_metric]], 
            color = .data[["model"]]
        )
    )
    if(!is.null(ass2d$benchmark) && !is.null(bm_data)){
        bm_alpha <- ifelse(
            ass2d$smooth_benchmark, 
            ass2d$alpha,
            1.
        )
        plt <- plt + ggplot2::geom_text(
            data = bm_data,
            mapping = ggplot2::aes(label = .data[["cutoff"]]),
            alpha = bm_alpha,
            size = ass2d$text_size,
            family = font_family
        )
    }
    plt <- plt +
        ggplot2::geom_point(alpha = ass2d$alpha) +
        ggplot2::labs(
            title = ass2d$title, 
            x = ass2d$x_lab, 
            y = ass2d$y_lab
        ) +
        ggplot2::scale_x_continuous(trans = ass2d$scale_x) + 
        ggplot2::scale_y_continuous(trans = ass2d$scale_y)
    if(!is.null(ass2d$hline))
        plt <- plt + do.call(ggplot2::geom_hline, ass2d$hline)
    if(!is.null(ass2d$vline))
        plt <- plt + do.call(ggplot2::geom_vline, ass2d$vline)
    if(!is.null(ass2d$text)){
        ass2d$text$family <- font_family
        plt <- plt + do.call(ggplot2::geom_label, ass2d$text)
    }
    if(!is.null(ass2d$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = ass2d$colors)
    }
    if(!is.null(ass2d$smooth_method)){
        plt <- plt + ggplot2::geom_smooth(
            method = ass2d$smooth_method,
            se = ass2d$smooth_se,
            formula = y ~ x
        )
        if(ass2d$smooth_benchmark && !is.null(bm_data)){
            plt <- plt + ggplot2::geom_smooth(
                data = bm_data,
                method = ass2d$smooth_method,
                se = ass2d$smooth_se,
                formula = y ~ x
            )
        }
    }
    if(!is.null(ass2d$theme)) plt <- plt + ass2d$theme

    if(ass2d$show_plots) print(plt)

    if(!is.null(ass2d$file)){
        if(!quiet)
            message(msg_prefix, "Saving 2D metric plot to ", ass2d$file)
        ggplot2::ggsave(
            ass2d$file, 
            plt, 
            width = ass2d$width, 
            height = ass2d$height, 
            units = ass2d$units,
            dpi = ass2d$dpi
        )

        # Save to csv (if wanted)
        if(ass2d$fellow_csv){
            csv_file <- stringr::str_replace(ass2d$file, "\\..+", ".csv")
            if(!quiet)
                message(msg_prefix, "Saving 2D metric table to ", csv_file)
            readr::write_csv(ass2d$data, csv_file)
        }
    }

    invisible(all_data)
}


#' @importFrom rlang .data
as2_plot_risk_scores <- function(self, private, data, model, quiet, msg_prefix){

    prep <- model$predict(
        data,
        quiet = quiet
    )
    pa <- intersect_by_names(prep[["predicted"]], prep[["actual"]], 
        rm_na = c(TRUE, TRUE))
    predicted <- pa[[1]]
    actual <- pa[[2]]

    tbl <- tibble::tibble(
        patient_id = names(predicted),
        rank = rank(-predicted),
        `risk score` = predicted,
        `true risk` = ifelse(unlist(actual) == 1, "high", "low")
    )

    plt <- ggplot2::ggplot(
        tbl,
        ggplot2::aes(
            x = .data[["rank"]], 
            y = .data[["risk score"]], 
            color = .data[["true risk"]]
            )
        ) +
        ggplot2::geom_point() +
        ggplot2::labs(
            title = self$title
        )
    if(!is.null(self$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = self$colors)
    }

    if(self$show_plots){
        print(plt)
    }

    if(!is.null(self$file)){
        if(!quiet)
            message(msg_prefix, "Saving scores plot to ", self$file)
        ggplot2::ggsave(self$file, plt, width = self$width, height = self$height, 
            units = self$units)

        if(self$fellow_csv){
            csv_file <- stringr::str_replace(self$file, "\\..+$", ".csv")
            if(!quiet)
                message(msg_prefix, "Saving scores table to ", csv_file)
            readr::write_csv(tbl, csv_file)
        }
    }
    invisible(tbl)
}