#' @title An R6 class to assess a model with one metric versus another.
#' @description Assess how well a model can predict time to event less than a certain
#' threshold with two metrics (like the ROC-AUC and the Brier score) to guide thresholding 
#' the continuous to binary model output.
#' @export
Ass2d <- R6::R6Class("Ass2d",
    public = list(
        
        #' @field x_metric,y_metric Metrics shown on both axes.
        x_metric = NULL,
        y_metric = NULL,
        #' @field pivot_time_cutoff Assess for predicting time to event less than 
        #' pivot_time_cutoff
        pivot_time_cutoff = NULL,
        #' @field benchmark Incorporate this benchmark (from pheno data) into 
        #' assessment.
        benchmark = NULL,
        #' @field file Store results in this file.
        file = NULL,
        #' @field lambda Assess model with this regularization parameter.
        lambda = NULL,
        #' @field ci_level For any confidence intervals calculated use this level.
        ci_level = NULL,
        #' @field fellow_csv Whether to also store the plotted data in a csv 
        #' including cutoffs.
        fellow_csv = NULL,
        #' @field scores_plot Whether to additionally plot rank versus model output score.
        scores_plot = NULL,
        #' @field show_plots Whether to print all generated plots.
        show_plots = NULL,
        #' @field title Plot title.
        title = NULL,
        #' @field x_lab,y_lab Axis labels.
        x_lab = NULL,
        y_lab = NULL,
        #' @field xlim,ylim Axis limits.
        xlim = NULL,
        ylim = NULL,
        #' @field smooth_method Method to smooth the plot.
        smooth_method = NULL,
        #' @field smooth_benchmark Whether to also smooth the benchmark.
        smooth_benchmark = NULL,
        #' @field smooth_se Whether to add standard error bands to the smoothed lines.
        smooth_se = NULL,
        #' @field scale_x,scale_y The scale of the axes.
        scale_x = NULL,
        scale_y = NULL,
        #' @field vline,hline Vertical/horizontal lines to be added to the plot.
        vline = NULL,
        hline = NULL,
        #' @field text_size The size of in-plot text.
        text_size = NULL,
        #' @field text Text label added to the plot.
        text = NULL,
        #' @field alpha The alpha value for the points in the plot.
        alpha = NULL,
        #' @field colors Colors used by ggplot2.
        colors = NULL,
        #' @field theme Theme applied to the plot.
        theme = NULL,
        #' @field width,height,units The size of the plot in units.
        width = NULL,
        height = NULL,
        units = NULL,
        #' @field dpi Plot resolution in dots per inch.
        dpi = NULL,

        #' @description Create a new Ass2d instance.
        #' @param x_metric string. The name of the performance measure to be plotted on the x-axis.
        #' All measures that can be passed to the `x.measure` parameter of [ROCR::performance()] are
        #' valid.
        #' @param y_metric string. The name of the performance measure to be plotted on the y-axis.
        #' 
        #' * All measures that can be passed to the `measure` parameter of [ROCR::performance()] are
        #' valid.
        #' * `"logrank"` for the p-values of a logrank test. In this case, `x_metric` must be
        #' `"prevalence"` or `"rpp"`.
        #' * `"precision_ci"` for the lower bound of the `ci_level` (see below) confidence interval
        #' of the precision according to Clopper-Pearson.
        #' @param pivot_time_cutoff numeric. Time-to-event threshold that divides samples into a
        #' high/low-risk (time to event below/above `pivot_time_cutoff`) group. Model performance
        #' will be measured in terms of how well it can separate these two groups as a binary
        #' classifier.
        #' @param lambda numeric or string. The regularization parameter of the model to
        #' evaluate.
        #' @param benchmark string. Incorporate this benchmark (from pheno data) into the 
        #' assessment if `benchmark` is not `NULL`.
        #' @param file string. If not `NULL`, store the resulting plot in this file.
        #' @param ci_level numeric in \[0, 1\]. The level used to calculate confidence intervals.
        #' @param fellow_csv logical. Whether to also store the plotted data in a csv including
        #' cutoffs.
        #' @param scores_plot logical. Display the ordered scores output by the model in a scatter 
        #' plot.
        #' @param show_plots logical. Whether to show the plots after creation in an interactive 
        #' graphics device.
        #' @param title string. The title of the plot. Default is `NULL`, i.e. no title.
        #' @param x_lab,y_lab string. Axis labels. Default is `x_metric` and `y_metric`, respectively.
        #' @param xlim,ylim numeric vector of length 2. The limits for both axes. Default is
        #' no contraints.
        #' @param smooth_method string or function. Smooth method to plot an additional smoothed graph.
        #' If `NULL` (the default), no smoothing. Else we pass `smooth_method` as the `method` 
        #' parameter to [`ggplot2::geom_smooth()`].
        #' @param smooth_benchmark logical. Whether to also smooth the benchmark data.
        #' @param smooth_se logical. Whether to add standard error bands to the smoothed lines.
        #' @param scale_x,scale_y string or transformation object (see [`scales::trans_new`] for the 
        #' latter). The scale of the axes, we will pass them to the `trans` paramter of 
        #' [`ggplot2::scale_x_continuous()`], [`ggplot2::scale_y_continuous()`], respectively. 
        #' @param vline,hline list or NULL. Vertical/horizontal lines to be added to the plot. A list
        #' holding the arguments to pass to [`ggplot2::geom_vline()`] and [`ggplot2::geom_hline()`],
        #' respectively. Default is `NULL`, no line.
        #' @param text_size numeric. The size of text passed to the `size` parameter of 
        #' `ggplot2::geom_text()`. The reason this attribute exists is you might want to have text 
        #' size different the global text size of the plot. See \code{vignette("ggplot2-specs", 
        #' package = "ggplot2")} for more details.
        #' @param text list or NULL. Text label added to the plot. A list holding the arguments to
        #' pass to [`ggplot2::geom_label()`]. Default is `NULL`, no text at all.
        #' @param alpha numeric in \[0, 1\]. The alpha value for the points and lines in the 
        #' plot.
        #' @param width numeric. The width of the plot in `units`.
        #' @param height numeric. The height of the plot in `units`.
        #' @param colors character vector. The colors to be used for the different models.
        #' Default is `NULL`, which means that the default colors of `ggplot2` will be used.
        #' @param theme S3 object inheriting from `"theme"` and `"gg"` (typically the return value of 
        #' [`ggplot2::theme()`] or a complete ggplot2 theme like [`ggplot2::theme_light()`]). 
        #' The theme of the plot. Default is `NULL`, which means the default theme of ggplot2.
        #' @param units string. The units of `width` and `height`. Default is `"in"` (inches).
        #' @param dpi numeric. Plot resolution in dots per inch.
        #' @return A new Ass2d object.
        initialize = function(
            x_metric,
            y_metric,
            pivot_time_cutoff,
            lambda = "lambda.min",
            benchmark = NULL,
            file = NULL,
            ci_level = 0.95,
            fellow_csv = FALSE,
            scores_plot = FALSE,
            show_plots = FALSE,
            title = NULL,
            x_lab = NULL,
            y_lab = NULL,
            xlim = c(-Inf, Inf),
            ylim = c(-Inf, Inf),
            smooth_method = NULL,
            smooth_benchmark = FALSE,
            smooth_se = FALSE,
            scale_x = "identity",
            scale_y = "identity",
            vline = NULL,
            hline = NULL,
            text_size = 2.5,
            text = NULL,
            alpha = 0.5,
            width = 7,
            height = 4,
            colors = NULL,
            theme = NULL,
            units = "in",
            dpi = 300
        ){
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
            stopifnot(is.character(file))
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
        },

        #' @description Assess a *single* model (with multiple splits) on a data set.
        #' @param data Data object. Assess on this data. Data must already be read in and 
        #' `cohort` attribute set.
        #' @param model ModelSpec object. Assess this model, multiple splits are supported.
        #' @param quiet logical. Whether to suppress messages.
        #' @param msg_prefix string. Prefix for messages. Default is `""`.
        #' @details We add the data underlying the plots to Ass2d object as a new attribute 
        #' named `data`.
        #' @export
        assess = function(
            data,
            model,
            quiet = FALSE,
            msg_prefix = ""
        ){
            directory <- dirname(self$file)
            if(!dir.exists(directory))
                dir.create(directory, recursive = TRUE)

            # Prepare, predict and calculate performance metric
            # (a) For model
            prep <- prepare_and_predict(
                expr_mat = data$expr_mat,
                pheno_tbl = data$pheno_tbl,
                data = data,
                model = model,
                lambda = self$lambda,
                pivot_time_cutoff = self$pivot_time_cutoff,
                benchmark_col = self$benchmark
            )
            calculate_2d_metric(
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
                ass2d = self,
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
        },

        #' @description Assess *multiple* models (with multiple splits) on a data set.
        #' @param data DataSpec object. Assess on this data set.
        #' @param model_list list of ModelSpec objects. Assess these models.
        #' We infer the `AssSpec2d` for the single plots in a reasonable way from it.
        #' @param model_tree_mirror character vector of length 2. Store the single plots in 
        #' the same directory as the corresponding model for the train cohort and the 
        #' comparison plot at `file`. When doing the analogous for the test data, replace 
        #' the first element of `model_tree_mirror` by the second element in every file path.
        #' @param comparison_plot logical. Whether to generate a plot holding the assessment 
        #' of all models. Default is `TRUE`.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @export
        assess_center = function(
            data,
            model_list,
            model_tree_mirror = c("models", "results"),
            comparison_plot = TRUE,
            quiet = FALSE
        ){
            cohorts <- data$cohort
            if(is.null(cohorts)) cohorts <- "test"

            perf_tbls <- list()

            data$read()
            expr_mat <- data$expr_mat
            pheno_tbl <- data$pheno_tbl

            if(!quiet) message("\nASSESSING ON ", data$name)
            for(cohort in cohorts){
                if(!quiet) message("# On ", cohort, " cohort")
                data$cohort <- cohort
                for(model in model_list){
                    if(!quiet) message("## ", model$name)
                    for(time_cutoff in model$time_cutoffs){
                        if(!quiet) message("### At time cutoff ", time_cutoff)
                        model_cutoff <- at_time_cutoff(model, time_cutoff)
                        this_as2 <- self$infer(
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
                        perf_tbls[[model_cutoff$name]] <- this_as2$data
                    }
                }
                cohort_as2 <- self$clone()
                cohort_as2$set("public", "data", dplyr::bind_rows(perf_tbls))
                if(comparison_plot){
                    if(cohort == "test")
                        cohort_as2$file <- mirror_path(
                            filepath = ass2d$file,
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
        }
    ),

    private = list(
        infer = function(
            model,
            data,
            model_tree_mirror
        ){
            # Prepare for assess_2d()
            this_as2 <- self$clone()
            this_as2$file <- file.path(
                model$directory,
                paste0(
                    self$x_metric, "_vs_", self$y_metric,
                    stringr::str_extract(self$file, "\\..+$")
                )
            )
            if(data$cohort == "test")
                this_as2$file <- mirror_path(
                    filepath = this_as2$file,
                    mirror = self$model_tree_mirror
                )
            this_as2$title <- paste0(
                data$name, " ", data$cohort, ", ",
                data$time_to_event_col, " < ", this_as2$pivot_time_cutoff
            )
            return(this_as2)
        }
    )      
)
