#' @importFrom rlang .data
plot_perf_metric <- function(
    perf_tbl,
    perf_plot_spec
){
    xlab <- "heya"
    plt <- ggplot2::ggplot(
        perf_tbl,
        ggplot2::aes(
            x = .data[[perf_plot_spec$x_metric]], 
            y = .data[[perf_plot_spec$y_metric]], 
            color = .data[["model"]]
            )
        ) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
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

    ggplot2::ggsave(
        perf_plot_spec$fname, 
        plt, 
        width = perf_plot_spec$width, 
        height = perf_plot_spec$height, 
        units = perf_plot_spec$units
    )

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

    return(tbl)
}


add_benchmark_perf_metric <- function(
    pheno_tbl,
    data_spec,
    perf_plot_spec,
    model_spec
){
    
    if(is.null(perf_plot_spec$benchmark))
        return(NULL)

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

    perf_tbl <- calculate_perf_metric(
        predicted = pred_act[[1]],
        actual = pred_act[[2]],
        perf_plot_spec = perf_plot_spec
    )

    return(perf_tbl)
}