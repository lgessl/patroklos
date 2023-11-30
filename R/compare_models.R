compare_models <- function(
    model_spec_list,
    data_spec_list,
    perf_plot_spec
){
    # If model names for legend are not present, automatically generate them 
    if(is.null(perf_plot_spec$model_names_in_legend)){
        perf_plot_spec$model_names_in_legend <- 
            sapply(model_spec_list, function(model_spec) model_spec$fit_fname)
    }
    perf_tbls <- list()
    if(length(data_spec_list) == 1L){
        data_spec_list <- rep(data_spec_list, length(model_spec_list))
    }

    # Collect tibbles, (optionally) individual plots, csv files
    for(i in seq_along(model_spec_list)){
        data_spec <- data_spec_list[[i]]
        model_spec <- model_spec_list[[i]]
        data <- read(data_spec)
        pred_y <- prepare_and_predict(
            expr_mat = data[["expr_mat"]],
            pheno_tbl = data[["pheno_tbl"]],
            data_spec = data_spec,
            model_spec = model_spec
        )
        perf_tbls[[i]] <- assess_model(
            predictions = pred_y[["predictions"]],
            y = pred_y[["y"]],
            perf_plot_spec = perf_plot_spec,
            model_spec = model_spec
        )
    }
    names(perf_tbls) <- perf_plot_spec$model_names_in_legend
    perf_tbl <- dplyr::bind_rows(perf_tbls, .id = "model")

    # Plot
    plt <- ggplot2::ggplot(
        perf_tbl,
        ggplot2::aes(
            x = .data[[perf_plot_spec$x_metric]], 
            y = .data[[perf_plot_spec$y_metric]], 
            color = model
            )
        ) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::ggtitle(perf_plot_spec$title) +
        ggplot2::xlab(perf_plot_spec$x_axis) +
        ggplot2::ylab(perf_plot_spec$y_axis)
    if(!is.null(perf_plot_spec$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = perf_plot_spec$colors)
    }
    if(perf_plot_spec$show_plots){
        print(plt)
    }
    ggplot2::ggsave(
        perf_plot_spec$fname, 
        plt, 
        width = perf_plot_spec$width, 
        height = perf_plot_spec$height, 
        units = perf_plot_spec$units
    )

    return(perf_tbl)
}

assess_model <- function(
    predictions,
    y,
    perf_plot_spec,
    model_spec
){
    # Extract most frequently used values
    x_metric <- perf_plot_spec$x_metric
    y_metric <- perf_plot_spec$y_metric
    plt_basename <- basename(perf_plot_spec$fname)

    # Calculate performance measures
    rocr_prediction <- ROCR::prediction(predictions, y)
    rocr_perf <- ROCR::performance(
        rocr_prediction,
        measure =  x_metric,
        x.measure = y_metric
    )
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
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]

    # Plot
    if(perf_plot_spec$also_single_plots){
        plt <- ggplot2::ggplot(tbl, ggplot2::aes(.data[[x_metric]], .data[[y_metric]])) +
            ggplot2::geom_line() +
            ggplot2::geom_point() +
            ggplot2::ggtitle(perf_plot_spec$title) +
            ggplot2::xlab(perf_plot_spec$x_axis) +
            ggplot2::ylab(perf_plot_spec$y_axis)
        if(!is.null(perf_plot_spec$colors)){
            plt <- plt + ggplot2::scale_color_manual(values = perf_plot_spec$colors)
        }

        plt_fname <- file.path(model_spec$save_dir, plt_basename)
        ggplot2::ggsave(
            plt_fname, 
            plt, 
            width = perf_plot_spec$width, 
            height = perf_plot_spec$height, 
            units = perf_plot_spec$units
        )
        if(perf_plot_spec$show_plots){
            print(plt)
        }
    }

    # store csv
    if(perf_plot_spec$single_csvs){
        csv_basename <- stringr::str_replace(plt_basename, "\\..+$", ".csv")
        csv_fname <- file.path(model_spec$save_dir, csv_basename)
        readr::write_csv(tbl, csv_fname)
    }

    return(tbl)
}