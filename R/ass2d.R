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
        #' @field data The data underlying the plots.
        data = NULL,

        #' @description Create a new Ass2d instance.
        #' @param x_metric string. The name of the performance measure to be plotted on the x-axis.
        #' All measures that can be passed to the `x.measure` parameter of [`ROCR::performance()`] are
        #' valid.
        #' @param y_metric string. The name of the performance measure to be plotted on the y-axis.
        #' 
        #' * All measures that can be passed to the `measure` parameter of [`ROCR::performance()`] are
        #' valid.
        #' * `"logrank"` for the p-values of a logrank test. In this case, `x_metric` must be
        #' `"prevalence"` or `"rpp"`.
        #' * `"precision_ci"` for the lower bound of the `ci_level` (see below) confidence interval
        #' of the precision according to Clopper-Pearson.
        #' @param pivot_time_cutoff numeric. Time-to-event threshold that divides samples into a
        #' high/low-risk (time to event below/above `pivot_time_cutoff`) group. Model performance
        #' will be measured in terms of how well it can separate these two groups as a binary
        #' classifier.
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
        )
            ass2d_initialize(self, private, x_metric, y_metric, pivot_time_cutoff, 
                benchmark, file, ci_level, fellow_csv, scores_plot, show_plots, title, 
                x_lab, y_lab, xlim, ylim, smooth_method, smooth_benchmark, smooth_se, 
                scale_x, scale_y, vline, hline, text_size, text, alpha, colors, theme, 
                width, height, units, dpi),

        #' @description Assess a *single* model (with multiple splits) on a data set.
        #' @param data Data object. Assess on this data. The `cohort` attribute 
        #' must be set.
        #' @param model Model object. Assess this model, multiple splits are supported.
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
        )
            ass2d_assess(self, private, data, model, quiet, msg_prefix),

        #' @description Plot rank versus model score on a data set. You can do this 
        #' on the fly, i.e. you actually initialize the Ass2d object for another 
        #' purpose and this method will infer reasonable defaults for the plot.
        #' @param data Data object. Plot the scores for this data.
        #' @param model Model object. Plot the scores for this model.
        #' @param quiet logical. Whether to suppress messages.
        #' @param msg_prefix string. Prefix for messages.
        plot_risk_scores = function(
            data,
            model,
            quiet = FALSE,
            msg_prefix = ""
        )
            as2_plot_risk_scores(self, private, data, model, quiet, msg_prefix),

        #' @description Assess *multiple* models (with multiple splits) on a data set.
        #' This is a wrapper around `assess()` and optionally `plot_risk_scores()` 
        #' for multiple models.
        #' @param data Data object. Assess on this data set.
        #' @param model_list list of Model objects. Assess these models.
        #' We infer the `AssSpec2d` for the single plots in a reasonable way from it.
        #' We skip models with the `continuous_output` attribute being FALSE (a plot 
        #' is an overkill in this case).
        #' @param model_tree_mirror character vector of length 2. Store the single plots in 
        #' the same directory as the corresponding model for the train cohort and the 
        #' comparison plot at `file`. When doing the analogous for the test data, replace 
        #' the first element of `model_tree_mirror` by the second element in every file path.
        #' @param risk_scores logical. Whether to also plot rank versus risk score 
        #' along with the assessment "on-the-fly".
        #' @param comparison_plot logical. Whether to generate a plot holding the assessment 
        #' of all models. Default is `TRUE`.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @export
        assess_center = function(
            data,
            model_list,
            model_tree_mirror = c("models", "results"),
            risk_scores = FALSE,
            comparison_plot = FALSE,
            quiet = FALSE
        )
            ass2d_assess_center(self, private, data, model_list, model_tree_mirror, 
                risk_scores, comparison_plot, quiet)
    ),

    private = list(
        # Part of assess method
        calculate_2d_metric = function(
            data,
            model,
            quiet = FALSE
        )
            ass2d_calculate_2d_metric(self, private, data, model, quiet),

        # Infer reasonable values for a new Ass2d object
        infer = function(
            model,
            data,
            model_tree_mirror
        )
            ass2d_infer(self, private, model, data, model_tree_mirror)
    )      
)
