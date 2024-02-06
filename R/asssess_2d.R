#' @title Assess a model on a single data set in terms how well it filters high-risk
#' samples in a 2d plot
#' @description Given a data set and a model, assess how well the model can filter
#' high-risk patients. This includes: 
#' 1. A 2D scatter plot of two performance metrics, namely the `x_metric` and `y_metric`
#' attributes of `ass_2d_spec`. 
#' 2. A plot of the risk scores of the model.
#' Supports multiple splits into train and test cohorts, but not multiple time cutoffs 
#' or multiple models. For the this, see [`assess_2d_center()`].
#' @param expr_mat numeric matrix. The expression matrix. Rows are samples, columns
#' are genes, their identifiers are given as row and column names, respectively.
#' @param pheno_tbl tibble. The phenotype table. Observations are samples, columns 
#' are variables.
#' @param data_spec DataSpec object specifying `expr_mat` and `pheno_tbl`. See the 
#' constructor `DataSpec()` for details.
#' @param model_spec ModelSpec object. The model to be assessed. See the constructor
#' `ModelSpec()` for details.
#' @param ass_2d_spec Ass2dSpec object. The specifications for the plot. See
#' the constructor `Ass2dSpec()` for details. The `pivot_time_cutoff` attribute of
#' `ass_2d_spec` will override the same in `model_spec` if both are given.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
#' @return A Ass2dSpec object. `ass_2d_spec` with one additional attribute named
#' `data` holding the data underlying the plots.
#' @details The Ass2dSpec class is tailored for this function, so see its constructor 
#' [`Ass2dSpec()`] for details.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
assess_2d <- function(
    expr_mat,
    pheno_tbl,
    data_spec,
    model_spec,
    ass_2d_spec,
    quiet = FALSE,
    msg_prefix = ""
){
    if(is.null(ass_2d_spec$pivot_time_cutoff)){
        ass_2d_spec$pivot_time_cutoff <- data_spec$pivot_time_cutoff
        if(is.null(ass_2d_spec$pivot_time_cutoff)){
            stop("You need to specify the `pivot_time_cutoff` atrribute in at least ",
                "one of `data_spec` and `ass_2d_spec`.")
        }
    }
    if(is.null(ass_2d_spec$title)){
        ass_2d_spec$title <- paste0(
            data_spec$name, " ", data_spec$cohort, ", time cutoff ", 
            ass_2d_spec$pivot_time_cutoff)
    }
    if(is.null(data_spec$cohort)){
        data_spec$cohort <- "test"
    }
    directory <- dirname(ass_2d_spec$fname)
    if(!dir.exists(directory))
        dir.create(directory, recursive = TRUE)

    # Prepare, predict and calculate performance metric
    # (a) For model
    prep <- prepare_and_predict(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data_spec = data_spec,
        model_spec = model_spec,
        lambda = ass_2d_spec$lambda,
        pivot_time_cutoff = ass_2d_spec$pivot_time_cutoff,
        benchmark_col = ass_2d_spec$benchmark
    )
    ass_2d_spec <- calculate_2d_metric(
        actual = prep[["actual"]],
        predicted = prep[["predicted"]],
        ass_2d_spec = ass_2d_spec,
        model_spec = model_spec,
        benchmark = prep[["benchmark"]],
        pheno_tbl = pheno_tbl,
        data_spec = data_spec
    )

    # Plot
    plot_2d_metric(
        ass_2d_spec = ass_2d_spec,
        quiet = quiet,
        msg_prefix = msg_prefix
    )

    if(ass_2d_spec$scores_plot){
        pps_scores <- ass_2d_spec
        pps_scores$title <- paste0(model_spec$name, " | ", ass_2d_spec$title)
        pps_scores$fname <- file.path(
            dirname(ass_2d_spec$fname),
            paste0("scores", stringr::str_extract(ass_2d_spec$fname, "\\..+$"))
        )
        plot_risk_scores(
            predicted = prep[["predicted"]],
            actual = prep[["actual"]],
            ass_2d_spec = pps_scores,
            ncol = model_spec$plot_ncols,
            msg_prefix = msg_prefix
        )
    }

    return(ass_2d_spec)
}