calculate_2d_metric <- function(
    actual,
    predicted,
    perf_plot_spec,
    model_spec,
    benchmark = NULL
){
    # Extract most frequently used values
    x_metric <- perf_plot_spec$x_metric
    y_metric <- perf_plot_spec$y_metric
    # Prepare for loop
    tbl_list <- list()
    aucs <- rep(NA, length(actual))
    estimate_list <- list(predicted, benchmark)
    if(is.null(benchmark)) estimate_list <- list(predicted)
    names(estimate_list) <- c(model_spec$name, perf_plot_spec$benchmark)

    for(i in model_spec$split_index){
        for(estimate_name in names(estimate_list)){
            estimate <- estimate_list[[estimate_name]]
            # Calculate performance measures
            rocr_prediction <- ROCR::prediction(estimate[[i]], actual[[i]])
            rocr_perf <- ROCR::performance(
                rocr_prediction,
                measure =  y_metric,
                x.measure = x_metric
            )
            if(estimate_name == model_spec$name){
                # Infer ROC-AUC
                aucs[i] <- ROCR::performance(
                    rocr_prediction,
                    measure = "auc"
                )@y.values[[1]]
            }
            # Store them in a tibble
            tbl <- tibble::tibble(
                rocr_perf@x.values[[1]],
                rocr_perf@y.values[[1]],
                rocr_perf@alpha.values[[1]]
            )
            names(tbl) <- c(
                x_metric, 
                y_metric,
                "cutoff"
            )
            tbl[["split"]] <- i
            tbl[["model"]] <- estimate_name
            tbl_list <- c(tbl_list, list(tbl))
        }
    }
    auc_avg <- mean(aucs, na.rm = TRUE)
    avg_string <- ifelse(length(aucs) > 1, "avg ", "")

    # Combine all splits to one tibble
    tbl <- dplyr::bind_rows(tbl_list)
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]
    perf_plot_spec$data <- tbl
    perf_plot_spec$title <- stringr::str_c(
        perf_plot_spec$title, ", ", 
        "ROC-AUC = ", round(auc_avg, 3)
    )

    return(perf_plot_spec)
}


#' @importFrom rlang .data
plot_2d_metric <- function(
    perf_plot_spec,
    quiet = FALSE
){
    # Extract
    data <- perf_plot_spec$data
    x_metric <- perf_plot_spec$x_metric
    y_metric <- perf_plot_spec$y_metric
    # Constraint data to range
    data <- data[
        data[[x_metric]] >= perf_plot_spec$xlim[1] &
        data[[x_metric]] <= perf_plot_spec$xlim[2],
    ]
    data <- data[
        data[[y_metric]] >= perf_plot_spec$ylim[1] &
        data[[y_metric]] <= perf_plot_spec$ylim[2],
    ]
    bm_data <- data[data[["model"]] == perf_plot_spec$benchmark, ]
    data <- data[data[["model"]] != perf_plot_spec$benchmark, ]

    plt <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes(
            x = .data[[x_metric]], 
            y = .data[[y_metric]], 
            color = .data[["model"]]
            )
        ) +
        ggplot2::geom_point(alpha = perf_plot_spec$alpha) +
        ggplot2::labs(
            title = perf_plot_spec$title, 
            x = perf_plot_spec$x_lab, 
            y = perf_plot_spec$y_lab
        )
    if(!is.null(perf_plot_spec$benchmark) && !is.null(bm_data)){
        bm_alpha <- ifelse(
            perf_plot_spec$smooth_benchmark, 
            perf_plot_spec$alpha,
            1.
        )
        plt <- plt + ggplot2::geom_point(
            data = bm_data,
            alpha = bm_alpha
        )
    }
    if(!is.null(perf_plot_spec$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = perf_plot_spec$colors)
    }
    if(!is.null(perf_plot_spec$smooth_method)){
        plt <- plt + ggplot2::geom_smooth(
            method = perf_plot_spec$smooth_method,
            se = FALSE,
            formula = y ~ x
        )
        if(perf_plot_spec$smooth_benchmark && !is.null(bm_data)){
            plt <- plt + ggplot2::geom_smooth(
                data = bm_data,
                method = perf_plot_spec$smooth_method,
                se = FALSE,
                formula = y ~ x
            )
        }
    }

    if(perf_plot_spec$show_plots){
        print(plt)
    }

    if(!quiet)
        message("Saving 2D metric plot to ", perf_plot_spec$fname)

    ggplot2::ggsave(
        perf_plot_spec$fname, 
        plt, 
        width = perf_plot_spec$width, 
        height = perf_plot_spec$height, 
        units = perf_plot_spec$units
    )

    # Save to csv (if wanted)
    if(perf_plot_spec$fellow_csv){
        csv_fname <- stringr::str_replace(perf_plot_spec$fname, "\\..+", ".csv")
        if(!quiet)
            message("Saving 2D metric table to ", csv_fname)
        readr::write_csv(perf_plot_spec$data, csv_fname)
    }

    return(plt)
}


#' @importFrom rlang .data
plot_risk_scores <- function(
    predicted,
    actual,
    perf_plot_spec,
    ncol = 2,
    quiet = FALSE
){
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
            ncol = ncol
        ) +
        ggplot2::geom_point() +
        ggplot2::labs(
            title = perf_plot_spec$title
        )
    if(!is.null(perf_plot_spec$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = perf_plot_spec$colors)
    }

    if(perf_plot_spec$show_plots){
        print(plt)
    }

    if(!quiet)
        message("Saving scores plot to ", perf_plot_spec$fname)
    ggplot2::ggsave(
        perf_plot_spec$fname, 
        plt, 
        width = perf_plot_spec$width, 
        height = perf_plot_spec$height, 
        units = perf_plot_spec$units
    )

    if(perf_plot_spec$fellow_csv){
        csv_fname <- stringr::str_replace(perf_plot_spec$fname, "\\..+", ".csv")
        if(!quiet)
            message("Saving scores table to ", csv_fname)
        readr::write_csv(tbl, csv_fname)
    }

    return(tbl)
}