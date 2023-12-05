assess_model <- function(
    data_spec,
    model_spec,
    perf_plot_spec
){
    # Read in
    data <- read(data_spec)
    expr_mat <- data[["expr_mat"]]
    pheno_tbl <- data[["pheno_tbl"]]

    # Prepare, predict and calcualte performance metric
    # (a) For model
    pred_act <- prepare_and_predict(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec
    )
    actual <- pred_act[["actual"]]
    predicted <- pred_act[["predicted"]]
    perf_tbl_model <- calculate_perf_metric(
        predicted = predicted,
        actual = actual,
        perf_plot_spec = perf_plot_spec
    )

    # (b) For benchmark (if given)
    perf_tbl_bm <- add_benchmark_perf_metric(
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        perf_plot_spec = perf_plot_spec,
        model_spec = model_spec
    )

    # Combine both into tibble
    perf_tbl_list <- list()
    perf_tbl_list[[model_spec$name]] <- perf_tbl_model
    if(!is.null(perf_tbl_bm)){
        perf_tbl_list[[perf_plot_spec$benchmark]] <- perf_tbl_bm
    }
    perf_tbl <- dplyr::bind_rows(perf_tbl_list, .id = "model")

    # Plot
    if(perf_plot_spec$also_single_plots){
        plot_perf_metric(
            perf_tbl = perf_tbl,
            perf_plot_spec = perf_plot_spec
        )
    }
    # Save to csv (if wanted)
    if(perf_plot_spec$single_csvs){
        csv_fname <- stringr::str_replace(perf_plot_spec$fname, "\\..+", ".csv")
        readr::write_csv(perf_tbl_model, csv_fname)
    }

    perf_tbls <- list("benchmark" = perf_tbl_bm, "model" = perf_tbl_model)
    return(perf_tbls)
}