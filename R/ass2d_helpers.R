#' @title Assess a model on a single data set in terms how well it filters high-risk
#' samples in a 2d plot
#' @description Given a data set and a model, assess how well the model can filter
#' high-risk patients. This includes: 
#' 1. A 2D scatter plot of two performance metrics, namely the `x_metric` and `y_metric`
#' attributes of `ass2d`. 
#' 2. A plot of the risk scores of the model.
#' Supports multiple splits into train and test cohorts, but not multiple time cutoffs 
#' or multiple models. For the this, see [`assess_2d_center()`].
#' @param expr_mat numeric matrix. The expression matrix. Rows are samples, columns
#' are genes, their identifiers are given as row and column names, respectively.
#' @param pheno_tbl tibble. The phenotype table. Observations are samples, columns 
#' are variables.
#' @param data DataSpec object specifying `expr_mat` and `pheno_tbl`. See the 
#' constructor `DataSpec()` for details.
#' @param model ModelSpec object. The model to be assessed. See the constructor
#' `ModelSpec()` for details.
#' @param ass2d AssSpec2d object. The specifications for the plot. See
#' the constructor `AssSpec2d()` for details. The `pivot_time_cutoff` attribute of
#' `ass2d` will override the same in `model` if both are given.
#' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
#' @return A AssSpec2d object. `ass2d` with one additional attribute named
#' `data` holding the data underlying the plots.
#' @details The AssSpec2d class is tailored for this function, so see its constructor 
#' [`AssSpec2d()`] for details.
#' @param msg_prefix string. Prefix for messages. Default is `""`.
#' @export
assess_2d <- function(
    expr_mat,
    pheno_tbl,
    data,
    model,
    ass2d,
    quiet = FALSE,
    msg_prefix = ""
){
    if(is.null(ass2d$pivot_time_cutoff)){
        ass2d$pivot_time_cutoff <- data$pivot_time_cutoff
        if(is.null(ass2d$pivot_time_cutoff)){
            stop("You need to specify the `pivot_time_cutoff` atrribute in at least ",
                "one of `data` and `ass2d`.")
        }
    }
    if(is.null(ass2d$title)){
        ass2d$title <- paste0(
            data$name, " ", data$cohort, ", time cutoff ", 
            ass2d$pivot_time_cutoff)
    }
    if(is.null(data$cohort)){
        data$cohort <- "test"
    }
    directory <- dirname(ass2d$file)
    if(!dir.exists(directory))
        dir.create(directory, recursive = TRUE)

    # Prepare, predict and calculate performance metric
    # (a) For model
    prep <- prepare_and_predict(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data = data,
        model = model,
        lambda = ass2d$lambda,
        pivot_time_cutoff = ass2d$pivot_time_cutoff,
        benchmark_col = ass2d$benchmark
    )
    ass2d <- calculate_2d_metric(
        actual = prep[["actual"]],
        predicted = prep[["predicted"]],
        ass2d = ass2d,
        model = model,
        benchmark = prep[["benchmark"]],
        pheno_tbl = pheno_tbl,
        data = data
    )

    # Plot
    plot_2d_metric(
        ass2d = ass2d,
        quiet = quiet,
        msg_prefix = msg_prefix
    )

    if(ass2d$scores_plot){
        pps_scores <- ass2d
        pps_scores$title <- paste0(model$name, " | ", ass2d$title)
        pps_scores$file <- file.path(
            dirname(ass2d$file),
            paste0("scores", stringr::str_extract(ass2d$file, "\\..+$"))
        )
        plot_risk_scores(
            predicted = prep[["predicted"]],
            actual = prep[["actual"]],
            ass2d = pps_scores,
            quiet = quiet,
            ncol = model$plot_ncols,
            msg_prefix = msg_prefix
        )
    }

    return(ass2d)
}

