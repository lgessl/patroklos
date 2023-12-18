#' @importFrom rlang .data
plot_perf_metric <- function(
    perf_plot_spec,
    quiet = FALSE
){
    perf_tbl <- perf_plot_spec$data
    plt <- ggplot2::ggplot(
        perf_tbl,
        ggplot2::aes(
            x = .data[[perf_plot_spec$x_metric]], 
            y = .data[[perf_plot_spec$y_metric]], 
            color = .data[["model"]]
            )
        ) +
        ggplot2::geom_line(alpha = perf_plot_spec$alpha) +
        ggplot2::geom_point(alpha = perf_plot_spec$alpha) +
        ggplot2::labs(
            title = perf_plot_spec$title, 
            x = perf_plot_spec$x_lab, 
            y = perf_plot_spec$y_lab
        )
    if(!is.null(perf_plot_spec$colors)){
        plt <- plt + ggplot2::scale_color_manual(values = perf_plot_spec$colors)
    }

    if(perf_plot_spec$show_plots){
        print(plt)
    }

    if(!quiet)
        message("Saving performance plot to ", perf_plot_spec$fname)
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
            message("Saving performance table to ", csv_fname)
        readr::write_csv(perf_tbl, csv_fname)
    }

    return(plt)
}


calculate_perf_metric <- function(
    predicted,
    actual,
    perf_plot_spec
){
    # Extract most frequently used values
    x_metric <- perf_plot_spec$x_metric
    y_metric <- perf_plot_spec$y_metric

    # Calculate performance measures
    rocr_prediction <- ROCR::prediction(predicted, actual)
    rocr_perf <- ROCR::performance(
        rocr_prediction,
        measure =  y_metric,
        x.measure = x_metric
    )
    # By default, also infer ROC-AUC
    roc_auc <- ROCR::performance(
        rocr_prediction,
        measure = "auc"
    )@y.values[[1]]

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

    # Guarantee availability
    any_na <- apply(tbl, 1, function(x) any(is.na(x)))
    tbl <- tbl[!any_na, ]

    perf_plot_spec$data <- tbl
    perf_plot_spec$title <- stringr::str_c(
        perf_plot_spec$title, 
        ", ROC-AUC = ", round(roc_auc, 3)
    )

    return(perf_plot_spec)
}


add_benchmark_perf_metric <- function(
    pheno_tbl,
    data_spec,
    perf_plot_spec,
    model_spec
){
    
    if(is.null(perf_plot_spec$benchmark))
        stop("Cannot calculate performance metric for benchmark classifier ",
            "if it is not provided as `benchmark` in `perf_plot_spec`.")

    bm_name <- perf_plot_spec$benchmark
    if(is.null(pheno_tbl[[bm_name]])){
        stop("There is no column named ", bm_name, " in ",
            file.path(data_spec$directory, data_spec$pheno_fname))
    }

    predicted <- pheno_tbl[[bm_name]]
    names(predicted) <- pheno_tbl[[data_spec$patient_id_col]]
    model_spec$response_type <- "binary" # always evaluate for discretized response
    actual <- generate_response(
        pheno_tbl = pheno_tbl, 
        data_spec = data_spec, 
        model_spec = model_spec
    )[, 1]

    pred_act <- intersect_by_names(
        predicted,
        actual,
        rm_na = TRUE
    )

    pps_bm <- calculate_perf_metric(
        predicted = pred_act[[1]],
        actual = pred_act[[2]],
        perf_plot_spec = perf_plot_spec
    )
    perf_plot_spec$bm_data <- pps_bm$data

    return(perf_plot_spec)
}


#' @importFrom rlang .data
plot_scores <- function(
    predicted,
    actual,
    perf_plot_spec,
    quiet = FALSE
){
    true_risk <- ifelse(actual == 1, "high", "low") |> as.factor()
    tbl <- tibble::tibble(
        patient_id = names(predicted),
        rank = rank(-predicted),
        `risk score` = predicted,
        `true risk` = true_risk
    )
    tbl <- tbl[order(tbl[["rank"]]), ]
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