ass2d_initialize <- function(self, private, x_metric, y_metric, pivot_time_cutoff, 
    lambda, benchmark, file, ci_level, fellow_csv, scores_plot, show_plots, title, 
    x_lab, y_lab, xlim, ylim, smooth_method, smooth_benchmark, smooth_se, 
    scale_x, scale_y, vline, hline, text_size, text, alpha, colors, theme, 
    width, height, units, dpi)
{
    if(is.null(x_lab)){
        x_lab <- x_metric
    }
    if(is.null(y_lab)){
        y_lab <- y_metric
    }
    if(y_metric == "logrank"){
        if(!any(stringr::str_detect(c("prevalence", "rpp"), x_metric)))
            stop("For `y_metric` = 'logrank', `x_metric` must be 'prevalence' or 'rpp'.")
    }
    stopifnot(is.character(file) || is.null(file))
    stopifnot(is.character(x_metric))
    stopifnot(is.character(y_metric))
    stopifnot(is.null(pivot_time_cutoff) || is.numeric(pivot_time_cutoff))
    stopifnot(is.character(lambda) || is.numeric(lambda))
    stopifnot(is.character(benchmark) || is.null(benchmark))
    stopifnot(is.numeric(ci_level) && ci_level >= 0 && ci_level <= 1)
    stopifnot(is.logical(fellow_csv))
    stopifnot(is.logical(scores_plot))
    stopifnot(is.logical(show_plots))
    stopifnot(is.character(title) || is.null(title))
    stopifnot(is.character(x_lab))
    stopifnot(is.character(y_lab))
    stopifnot(is.numeric(xlim) || is.null(xlim) || 
        xlim[1] >= xlim[2] || ylim[1] >= ylim[2])
    stopifnot(is.numeric(ylim) || is.null(ylim))
    stopifnot(is.character(smooth_method) || is.null(smooth_method) ||
        is.function(smooth_method))
    stopifnot(is.logical(smooth_benchmark))
    stopifnot(is.logical(smooth_se))
    stopifnot(is.null(hline) || is.list(hline))
    stopifnot(is.null(vline) || is.list(vline))
    stopifnot(is.numeric(text_size) && text_size > 0)
    stopifnot(is.null(text) || is.list(text))
    stopifnot(is.numeric(alpha) && alpha >= 0 && alpha <= 1)
    stopifnot(is.character(colors) || is.null(colors))
    stopifnot(is.numeric(width) && width > 0)
    stopifnot(is.numeric(height) && height > 0)
    stopifnot((inherits(theme, "theme") && inherits(theme, "gg")) || is.null(theme))
    stopifnot(is.character(units))
    stopifnot(is.numeric(dpi) && dpi > 0)

    self$x_metric <- x_metric
    self$y_metric <- y_metric
    self$pivot_time_cutoff <- pivot_time_cutoff
    self$lambda <- lambda
    self$benchmark <- benchmark
    self$file <- file
    self$ci_level <- ci_level
    self$fellow_csv <- fellow_csv
    self$scores_plot <- scores_plot
    self$show_plots <- show_plots
    self$title <- title
    self$x_lab <- x_lab
    self$y_lab <- y_lab
    self$xlim <- xlim
    self$ylim <- ylim
    self$smooth_method <- smooth_method
    self$smooth_benchmark <- smooth_benchmark
    self$smooth_se <- smooth_se
    self$scale_x <- scale_x
    self$scale_y <- scale_y
    self$vline <- vline
    self$hline <- hline
    self$text_size <- text_size
    self$text <- text
    self$alpha <- alpha
    self$colors <- colors
    self$theme <- theme
    self$width <- width
    self$height <- height
    self$units <- units
    self$dpi <- dpi

    invisible(self)
}

ass2d_assess <- function(self, private, data, model, quiet, msg_prefix){

    directory <- dirname(self$file)
    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    # Has at_time_cutoff() been called already? Skip in this case
    if(!dir.exists(directory))
        dir.create(directory, recursive = TRUE)
    private$calculate_2d_metric(
        data = data, 
        model = model
    )
    plot_2d_metric(
        ass2d = self,
        quiet = quiet,
        msg_prefix = msg_prefix
    )

    invisible(self)
}

ass2d_assess_center <- function(self, private, data, model_list, model_tree_mirror, 
    risk_scores, comparison_plot, quiet){
    
    cohorts <- data$cohort
    if(is.null(cohorts)) cohorts <- "test"
    n_assess <- sum(sapply(model_list, function(m) length(m$time_cutoffs)))
    perf_tbls <- vector("list", length(cohorts)*n_assess)
    data$read()

    if(!quiet) message("\nASSESSING ON ", data$name)
    i <- 1
    for(cohort in cohorts){
        if(!quiet) message("# On ", cohort, " cohort")
        data$cohort <- cohort
        for(model in model_list){
            if(!quiet) message("## ", model$name)
            for(time_cutoff in model$time_cutoffs){
                if(!quiet) message("### At time cutoff ", time_cutoff)
                model_cutoff <- model$at_time_cutoff(time_cutoff)
                this_as2 <- private$infer(
                    model = model_cutoff,
                    data = data,
                    model_tree_mirror = model_tree_mirror
                )
                this_as2$assess(
                    data = data,
                    model = model_cutoff,
                    quiet = quiet,
                    msg_prefix = "#### "
                )
                perf_tbls[[i]] <- this_as2$data
                if(risk_scores)
                    this_as2$plot_risk_scores(
                        data = data,
                        model = model_cutoff,
                        quiet = quiet,
                        msg_prefix = "#### "
                    )
                i <- i+1
            }
        }
        cohort_as2 <- self$clone()
        cohort_as2$data <- dplyr::bind_rows(perf_tbls)
        if(comparison_plot){
            if(cohort == "test")
                cohort_as2$file <- mirror_path(
                    filepath = self$file,
                    mirror = model_tree_mirror
                )
            if(is.null(cohort_as2$title))
                cohort_as2$title <- paste0(
                    data$name, " ", data$cohort, ", ", data$time_to_event_col,
                    " < ", cohort_as2$pivot_time_cutoff
                )
            plot_2d_metric(
                ass2d = cohort_as2,
                quiet = TRUE
            )
            if(!quiet)
                message("# Saving comparative performance plot to ", cohort_as2$file)
        }
    }
    data$cohort <- cohorts # No side effects
}

ass2d_infer <- function(self, private, model, data, model_tree_mirror){
    this_as2 <- self$clone()
    file_ext <- ".jpeg"
    if(!is.null(self$file))
        file_ext <- stringr::str_extract(self$file, "\\..+$")
    this_as2$file <- file.path(
        model$directory,
        paste0(self$x_metric, "_vs_", self$y_metric, file_ext)
    )
    if(data$cohort == "test")
        this_as2$file <- mirror_path(
            filepath = this_as2$file,
            mirror = model_tree_mirror
        )
    this_as2$title <- paste0(
        data$name, " ", data$cohort, ", ",
        data$time_to_event_col, " < ", this_as2$pivot_time_cutoff
    )
    return(this_as2)
}