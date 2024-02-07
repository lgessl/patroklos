#' @importFrom rlang .data
plot_2d_metric <- function(
    ass_2d_spec,
    quiet = FALSE,
    msg_prefix = ""
){
    # Extract
    data <- ass_2d_spec$data |> dplyr::distinct()
    x_metric <- ass_2d_spec$x_metric
    y_metric <- ass_2d_spec$y_metric
    # Constraint data to range
    data <- data[
        data[[x_metric]] >= ass_2d_spec$xlim[1] &
        data[[x_metric]] <= ass_2d_spec$xlim[2],
    ]
    data <- data[
        data[[y_metric]] >= ass_2d_spec$ylim[1] &
        data[[y_metric]] <= ass_2d_spec$ylim[2],
    ]
    bm_data <- NULL
    if(!is.null(ass_2d_spec$benchmark)){
        bm_data <- data[data[["model"]] == ass_2d_spec$benchmark, ]
        data <- data[data[["model"]] != ass_2d_spec$benchmark, ]
    }

    plt <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = .data[[x_metric]], 
            y = .data[[y_metric]], 
            color = .data[["model"]]
            )
        ) +
        ggplot2::geom_point(alpha = ass_2d_spec$alpha) +
        ggplot2::labs(
            title = ass_2d_spec$title, 
            x = ass_2d_spec$x_lab, 
            y = ass_2d_spec$y_lab
        ) +
        ggplot2::scale_x_continuous(trans = ass_2d_spec$scale_x) + 
        ggplot2::scale_y_continuous(trans = ass_2d_spec$scale_y)
    if(!is.null(ass_2d_spec$hline))
        plt <- plt + do.call(ggplot2::geom_hline, ass_2d_spec$hline)
    if(!is.null(ass_2d_spec$vline))
        plt <- plt + do.call(ggplot2::geom_vline, ass_2d_spec$vline)
    if(!is.null(ass_2d_spec$text))
        plt <- plt + do.call(ggplot2::geom_label, ass_2d_spec$text)
    if(!is.null(ass_2d_spec$benchmark) && !is.null(bm_data)){
        bm_alpha <- ifelse(
            ass_2d_spec$smooth_benchmark, 
            ass_2d_spec$alpha,
            1.
        )
        plt <- plt + ggplot2::geom_point(
            data = bm_data,
            alpha = bm_alpha
        )
    }
    if(!is.null(ass_2d_spec$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = ass_2d_spec$colors)
    }
    if(!is.null(ass_2d_spec$smooth_method)){
        plt <- plt + ggplot2::geom_smooth(
            method = ass_2d_spec$smooth_method,
            se = FALSE,
            formula = y ~ x
        )
        if(ass_2d_spec$smooth_benchmark && !is.null(bm_data)){
            plt <- plt + ggplot2::geom_smooth(
                data = bm_data,
                method = ass_2d_spec$smooth_method,
                se = FALSE,
                formula = y ~ x
            )
        }
    }

    if(ass_2d_spec$show_plots){
        print(plt)
    }

    if(!quiet)
        message(msg_prefix, "Saving 2D metric plot to ", ass_2d_spec$file)

    ggplot2::ggsave(
        ass_2d_spec$file, 
        plt, 
        width = ass_2d_spec$width, 
        height = ass_2d_spec$height, 
        units = ass_2d_spec$units
    )

    # Save to csv (if wanted)
    if(ass_2d_spec$fellow_csv){
        csv_file <- stringr::str_replace(ass_2d_spec$file, "\\..+", ".csv")
        if(!quiet)
            message(msg_prefix, "Saving 2D metric table to ", csv_file)
        readr::write_csv(ass_2d_spec$data, csv_file)
    }

    return(plt)
}


#' @importFrom rlang .data
plot_risk_scores <- function(
    predicted,
    actual,
    ass_2d_spec,
    ncol = 2,
    quiet = FALSE,
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
    ass_2d_spec$height <- nrow * ass_2d_spec$height
    ass_2d_spec$width <- ncol * ass_2d_spec$width

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
            title = ass_2d_spec$title
        )
    if(!is.null(ass_2d_spec$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = ass_2d_spec$colors)
    }

    if(ass_2d_spec$show_plots){
        print(plt)
    }

    if(!quiet)
        message(msg_prefix, "Saving scores plot to ", ass_2d_spec$file)
    ggplot2::ggsave(
        ass_2d_spec$file, 
        plt, 
        width = ass_2d_spec$width, 
        height = ass_2d_spec$height, 
        units = ass_2d_spec$units
    )

    if(ass_2d_spec$fellow_csv){
        csv_file <- stringr::str_replace(ass_2d_spec$file, "\\..+$", ".csv")
        if(!quiet)
            message(msg_prefix, "Saving scores table to ", csv_file)
        readr::write_csv(tbl, csv_file)
    }

    return(tbl)
}