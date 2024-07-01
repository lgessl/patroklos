ass2d_initialize <- function(self, private, x_metric, y_metric, 
    benchmark, file, ci_level, fellow_csv, show_plots, title, 
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
    stopifnot(is.character(benchmark) || is.null(benchmark))
    stopifnot(is.numeric(ci_level) && ci_level >= 0 && ci_level <= 1)
    stopifnot(is.logical(fellow_csv))
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
    self$benchmark <- benchmark
    self$file <- file
    self$ci_level <- ci_level
    self$fellow_csv <- fellow_csv
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

    if (!is.null(self$file)) {
        directory <- dirname(self$file)
        if(!dir.exists(directory))
            dir.create(directory, recursive = TRUE)
    }
    if(is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()
    if (self$x_metric == "rank" && self$y_metric == "risk score")
        private$plot_risk_scores(data, model, quiet, msg_prefix)
    else {
        private$calculate_2d_metric(
            data = data, 
            model = model,
            quiet = quiet
        )
        plot_2d_metric(
            ass2d = self,
            quiet = quiet,
            msg_prefix = msg_prefix
        )
    }
}

ass2d_assess_center <- function(self, private, data, model_list, model_tree_mirror, 
    comparison_plot, quiet){
    
    if (self$x_metric == "rank" && self$y_metric == "risk score") {
        comparison_plot <- FALSE
        if (!quiet) message("Setting comparison plot to FALSE for ", 
            "rank versus risk score.")
    }
    if (is.null(data$expr_mat) || is.null(data$pheno_tbl)) data$read()

    if(!quiet) message("Assessment center on ", data$name, " open")
    perf_tbls <- vector("list", length(model_list))

    # Single plots
    for (i in seq_along(model_list)) {
        model <- model_list[[i]]
        if (!is.null(model$continuous_output) && !model$continuous_output)
            next
        if(!quiet) message("** ", model$name)
        this_as2 <- private$infer(
            model = model,
            data = data,
            model_tree_mirror = model_tree_mirror
        )
        perf_tbls[[i]] <- this_as2$assess(
            data = data,
            model = model,
            quiet = quiet,
            msg_prefix = "**** "
        )
    }
    perf_tbls <- perf_tbls[!sapply(perf_tbls, is.null)]

    # Comparison plot
    self$data <- dplyr::bind_rows(perf_tbls)
    if(comparison_plot){
        if (data$cohort == "test")
            self$file <- stringr::str_replace(self$file, model_tree_mirror[1],
                model_tree_mirror[2])
        if(is.null(self$title))
            self$title <- paste0(
                data$name, " ", data$cohort, ", ", data$time_to_event_col,
                " < ", data$pivot_time_cutoff
            )
        plot_2d_metric(
            ass2d = self,
            quiet = TRUE
        )
        if(!quiet)
            message("* Saving comparison plot to ", self$file)
    }
    perf_tbl <- self$data
    self$data <- NULL
    return(perf_tbl)
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
    this_as2$file <- stringr::str_replace_all(this_as2$file, "\\s+", "_")
    if (!is.null(model_tree_mirror) && data$cohort == "test")
        this_as2$file <- stringr::str_replace(this_as2$file, model_tree_mirror[1], 
            model_tree_mirror[2])
    this_as2$title <- paste0(
        data$name, " ", data$cohort, ", ",
        data$time_to_event_col, " < ", data$pivot_time_cutoff
    )
    return(this_as2)
}