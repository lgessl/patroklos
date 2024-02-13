#' @importFrom rlang .data
plot_2d_metric <- function(
    ass_spec_2d,
    quiet = FALSE,
    msg_prefix = ""
){
    # Extract
    data <- ass_spec_2d$data |> dplyr::distinct()
    # For legend order: bechmark last
    if(!is.null(ass_spec_2d$benchmark)){
        model_levels <- unique(data[["model"]])
        model_levels <- sort(model_levels[model_levels != ass_spec_2d$benchmark])
        model_levels <- c(ass_spec_2d$benchmark, model_levels)
        data[["model"]] <- factor(data[["model"]], levels = model_levels)
    }
    x_metric <- ass_spec_2d$x_metric
    y_metric <- ass_spec_2d$y_metric
    # Constraint data to range
    full_data <- data
    data <- data[
        data[[x_metric]] >= ass_spec_2d$xlim[1] &
        data[[x_metric]] <= ass_spec_2d$xlim[2],
    ]
    data <- data[
        data[[y_metric]] >= ass_spec_2d$ylim[1] &
        data[[y_metric]] <= ass_spec_2d$ylim[2],
    ]
    bm_data <- NULL
    if(!is.null(ass_spec_2d$benchmark)){
        bm_data <- data[data[["model"]] == ass_spec_2d$benchmark, ]
        data <- data[data[["model"]] != ass_spec_2d$benchmark, ]
    }
    plt <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = .data[[x_metric]], 
            y = .data[[y_metric]], 
            color = .data[["model"]]
        )
    )
    if(!is.null(ass_spec_2d$benchmark) && !is.null(bm_data)){
        bm_alpha <- ifelse(
            ass_spec_2d$smooth_benchmark, 
            ass_spec_2d$alpha,
            1.
        )
        plt <- plt + ggplot2::geom_text(
            data = bm_data,
            mapping = ggplot2::aes(label = .data[["cutoff"]]),
            alpha = bm_alpha,
            size = ass_spec_2d$text_size
        )
    }
    plt <- plt +
        ggplot2::geom_point(alpha = ass_spec_2d$alpha) +
        ggplot2::labs(
            title = ass_spec_2d$title, 
            x = ass_spec_2d$x_lab, 
            y = ass_spec_2d$y_lab
        ) +
        ggplot2::scale_x_continuous(trans = ass_spec_2d$scale_x) + 
        ggplot2::scale_y_continuous(trans = ass_spec_2d$scale_y)
    if(!is.null(ass_spec_2d$hline))
        plt <- plt + do.call(ggplot2::geom_hline, ass_spec_2d$hline)
    if(!is.null(ass_spec_2d$vline))
        plt <- plt + do.call(ggplot2::geom_vline, ass_spec_2d$vline)
    if(!is.null(ass_spec_2d[["text"]]))
        plt <- plt + do.call(ggplot2::geom_label, ass_spec_2d[["text"]])
    if(!is.null(ass_spec_2d$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = ass_spec_2d$colors)
    }
    if(!is.null(ass_spec_2d$smooth_method)){
        plt <- plt + ggplot2::geom_smooth(
            method = ass_spec_2d$smooth_method,
            se = ass_spec_2d$smooth_se,
            formula = y ~ x
        )
        if(ass_spec_2d$smooth_benchmark && !is.null(bm_data)){
            plt <- plt + ggplot2::geom_smooth(
                data = bm_data,
                method = ass_spec_2d$smooth_method,
                se = ass_spec_2d$smooth_se,
                formula = y ~ x
            )
        }
    }

    if(ass_spec_2d$show_plots){
        print(plt)
    }

    if(!quiet)
        message(msg_prefix, "Saving 2D metric plot to ", ass_spec_2d$file)

    ggplot2::ggsave(
        ass_spec_2d$file, 
        plt, 
        width = ass_spec_2d$width, 
        height = ass_spec_2d$height, 
        units = ass_spec_2d$units
    )

    # Save to csv (if wanted)
    if(ass_spec_2d$fellow_csv){
        csv_file <- stringr::str_replace(ass_spec_2d$file, "\\..+", ".csv")
        if(!quiet)
            message(msg_prefix, "Saving 2D metric table to ", csv_file)
        readr::write_csv(ass_spec_2d$data, csv_file)
    }

    return(plt)
}


#' @importFrom rlang .data
plot_risk_scores <- function(
    predicted,
    actual,
    ass_spec_2d,
    quiet,
    ncol = 2,
    msg_prefix = ""
){
    # Get rid of NAs
    for(i in seq_along(predicted)){
        pred_act <- intersect_by_names(
            predicted[[i]], 
            actual[[i]], 
            rm_na = TRUE
        )
        predicted[[i]] <- pred_act[[1]]
        actual[[i]] <- pred_act[[2]]
    }
    tbl <- tibble::tibble(
        patient_id = names(unlist(predicted)),
        rank = unlist(lapply(predicted, function(x) rank(-x))),
        split = rep(
            which(!sapply(predicted, is.null)), 
            sapply(predicted, length)
        ),
        `risk score` = unlist(predicted),
        `true risk` = ifelse(unlist(actual) == 1, "high", "low")
    )
    tbl[["split"]] <- paste0("Split ", tbl[["split"]])

    n_split <- sum(!sapply(predicted, is.null))
    nrow <- ceiling(n_split/ncol)
    ass_spec_2d$height <- nrow * ass_spec_2d$height
    ass_spec_2d$width <- ncol * ass_spec_2d$width

    plt <- ggplot2::ggplot(
        tbl,
        ggplot2::aes(
            x = .data[["rank"]], 
            y = .data[["risk score"]], 
            color = .data[["true risk"]]
            )
        ) +
        ggplot2::facet_wrap(
            facets = ggplot2::vars(.data[["split"]]),
            ncol = ncol,
            scales = "free"
        ) +
        ggplot2::geom_point() +
        ggplot2::labs(
            title = ass_spec_2d$title
        )
    if(!is.null(ass_spec_2d$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = ass_spec_2d$colors)
    }

    if(ass_spec_2d$show_plots){
        print(plt)
    }

    if(!quiet)
        message(msg_prefix, "Saving scores plot to ", ass_spec_2d$file)
    ggplot2::ggsave(
        ass_spec_2d$file, 
        plt, 
        width = ass_spec_2d$width, 
        height = ass_spec_2d$height, 
        units = ass_spec_2d$units
    )

    if(ass_spec_2d$fellow_csv){
        csv_file <- stringr::str_replace(ass_spec_2d$file, "\\..+$", ".csv")
        if(!quiet)
            message(msg_prefix, "Saving scores table to ", csv_file)
        readr::write_csv(tbl, csv_file)
    }

    return(tbl)
}