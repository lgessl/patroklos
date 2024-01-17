#' @title Assess a model on a single data set in terms how well it filters 
#' high-risk patients
#' @description Given a data set and a model, assess how well the model can filter
#' high-risk patients. This includes: 
#' 1. A 2D scatter plot of two performance metrics, namely the `x_metric` and `y_metric`
#' attributes of `perf_plot_spec`. 
#' 2. A plot of the risk scores of the model.
#' Supports multiple splits into train and test cohorts, but not multiple time cutoffs 
#' or multiple models. For the this, see [`assessment_center()`].
#' @param expr_mat numeric matrix. The expression matrix. Rows are samples, columns
#' are genes, their identifiers are given as row and column names, respectively.
#' @param pheno_tbl tibble. The phenotype table. Observations are samples, columns 
#' are variables.
#' @param data_spec DataSpec object specifying `expr_mat` and `pheno_tbl`. See the 
#' constructor `DataSpec()` for details.
#' @param model_spec ModelSpec object. The model to be assessed. See the constructor
#' `ModelSpec()` for details.
#' @param perf_plot_spec PerfPlotSpec object. The specifications for the plot. See
#' the constructor `PerfPlotSpec()` for details. The `pivot_time_cutoff` attribute of
#' `perf_plot_spec` will override the same in `model_spec` if both are given.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
#' @return A PerfPlotSpec object. `perf_plot_spec` with one additional attribute named
#' `data` holding the data underlying the plots.
#' @details The PerfPlotSpec class is tailored for this function, so see its constructor 
#' [`PerfPlotSpec()`] for details.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
assess_model <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec,
    perf_plot_spec,
    quiet = FALSE,
    msg_prefix = ""
){
    if(is.null(perf_plot_spec$pivot_time_cutoff)){
        perf_plot_spec$pivot_time_cutoff <- data_spec$pivot_time_cutoff
        if(is.null(perf_plot_spec$pivot_time_cutoff)){
            stop("You need to specify the `pivot_time_cutoff` atrribute in at least ",
                "one of `data_spec` and `perf_plot_spec`.")
        }
    }
    if(is.null(perf_plot_spec$title)){
        perf_plot_spec$title <- paste0(
            data_spec$name, " ", data_spec$cohort, ", time cutoff ", 
            perf_plot_spec$pivot_time_cutoff)
    }
    if(is.null(data_spec$cohort)){
        data_spec$cohort <- "test"
    }
    directory <- dirname(perf_plot_spec$fname)
    if(!dir.exists(directory))
        dir.create(directory, recursive = TRUE)

    # Prepare, predict and calculate performance metric
    # (a) For model
    prep <- prepare_and_predict(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec,
        lambda = perf_plot_spec$lambda,
        perf_plot_spec = perf_plot_spec
    )
    perf_plot_spec <- calculate_2d_metric(
        actual = prep[["actual"]],
        predicted = prep[["predicted"]],
        perf_plot_spec = perf_plot_spec,
        model_spec = model_spec,
        benchmark = prep[["benchmark"]]
    )

    # Plot
    plot_2d_metric(
        perf_plot_spec = perf_plot_spec,
        quiet = quiet,
        msg_prefix = msg_prefix
    )

    if(perf_plot_spec$scores_plot){
        pps_scores <- perf_plot_spec
        pps_scores$title <- paste0(model_spec$name, " | ", perf_plot_spec$title)
        pps_scores$fname <- file.path(
            dirname(perf_plot_spec$fname),
            paste0("scores", stringr::str_extract(perf_plot_spec$fname, "\\..+$"))
        )
        plot_risk_scores(
            predicted = prep[["predicted"]],
            actual = prep[["actual"]],
            perf_plot_spec = pps_scores,
            msg_prefix = msg_prefix
        )
    }

    return(perf_plot_spec)
}