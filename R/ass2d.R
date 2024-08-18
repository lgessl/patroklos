#' @title An R6 class to assess a model with continuous output with one metric versus another.
#' @description Assess how well a model can predict time to event less than a certain
#' threshold with two metrics by thresholding the continuous output in every possible way and 
#' then plot one metric for a binary classifier against another one.
#' @export
Ass2d <- R6::R6Class("Ass2d",
    public = list(
        
        #' @field x_metric,y_metric Metrics shown on both axes.
        x_metric = NULL,
        y_metric = NULL,
        #' @field confidence_level Confidence level gamma, e.g., for confidence intervals.
        confidence_level = NULL,
        #' @field x_lab,y_lab Axis labels.
        x_lab = NULL,
        y_lab = NULL,
        #' @field xlim,ylim Axis limits.
        xlim = NULL,
        ylim = NULL,
        #' @field scale_x,scale_y The scale of the axes.
        scale_x = NULL,
        scale_y = NULL,
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

        #' @description Construct an `Ass2d` instance.
        #' @param x_metric,y_metric string. Metric shown on the x and y axis. 
        #' Regarding the choices: 
        #' * Every combination of (`measure`, `x.measure`) one can pass to 
        #'  [`ROCR::performance()`].
        #' * If `y_metric` is `"logrank"` or `"precision_ci"` (the lower limit of 
        #'  the `confidence_level` confidence interval of the precision, see below), 
        #' `x_metric` must be `"prevalence"` or `"rpp"` (rate of positive 
        #' predictions).
        #' * For (`x_metric`, `y_metric`) = (`"rank"`, "`risk score`"), plot the 
        #' rank of every risk score against the risk score and color by true 
        #' risk.
        #' @param confidence_level numeric in \[0, 1\]. Confidence level gamma for confidence 
        #' intervals.
        #' @param x_lab,y_lab string. Axis labels. Default is `x_metric` and `y_metric`. 
        #' @param xlim,ylim numeric vector of length 2. The limits for both axes. Default is
        #' no limits
        #' @param scale_x,scale_y string or transformation object (see [`scales::trans_new`] for the 
        #' latter). The scale of the axes, we will pass them to the `trans` parameter of 
        #' [`ggplot2::scale_x_continuous()`], [`ggplot2::scale_y_continuous()`], respectively. 
        #' @param alpha numeric in \[0, 1\]. The alpha value for the points and lines in the 
        #' plot.
        #' @param width numeric. The width of the plot in `units`.
        #' @param height numeric. The height of the plot in `units`.
        #' @param colors character vector. The colors to be used for the different models.
        #' Default is `NULL`, which means that the default colors of `ggplot2` will be used.
        #' @param theme S3 object inheriting from `"theme"` and `"gg"` (typically the return value 
        #' of [`ggplot2::theme()`] or a complete ggplot2 theme like [`ggplot2::theme_light()`]). 
        #' The theme of the plot. Default is `ggplot2::theme_minimal()`.
        #' @param units string. The units of `width` and `height`. Default is `"in"` (inches).
        #' @param dpi numeric. Plot resolution in dots per inch.
        #' @return A new `Ass2d` object.
        initialize = function(
            x_metric,
            y_metric,
            confidence_level = 0.95,
            x_lab = NULL,
            y_lab = NULL,
            xlim = c(-Inf, Inf),
            ylim = c(-Inf, Inf),
            scale_x = "identity",
            scale_y = "identity",
            alpha = 1,
            colors = NULL,
            theme = ggplot2::theme_minimal(),
            width = 7,
            height = 4,
            units = "in",
            dpi = 300
        )
            ass2d_initialize(self, private, x_metric, y_metric, confidence_level, 
                x_lab, y_lab, xlim, ylim, scale_x, scale_y, alpha, 
                colors, theme, width, height, units, dpi),

        #' @description Assess a *single* model.
        #' @param data Data object. Assess on this data. The `cohort` attribute 
        #' must be set.
        #' @param model Model object. Assess this model.
        #' @param return_type string. Either "ggplot" or "tibble". See return section for details.
        #' @param file string. If not `NULL` and `return_type == "ggplot`, store the resulting plot 
        #' in this file.
        #' @param fellow_csv logical. If `TRUE`, `file` is not `NULL` and `return_type == "ggplot"`,
        #' store the plotted data in a csv named `file` with replaced file extension.
        #' @param quiet logical. Whether to suppress messages.
        #' @param msg_prefix string. Prefix for messages. Default is `""`.
        #' @return ggplot object if `return_type == "ggplot"` or the tibble underlying the plot if 
        #' `return_type == "tibble"`.
        #' @export
        assess = function(
            data,
            model,
            return_type = "ggplot",
            file = NULL,
            fellow_csv = FALSE,
            quiet = FALSE,
            msg_prefix = ""
        )
            ass2d_assess(self, private, data, model, return_type, file, fellow_csv, quiet, 
                msg_prefix),

        #' @description Wrap `assess()` to assess *multiple* models and store the result.
        #' @param data Data object. Assess on this data. The `cohort` attribute must be set.
        #' @param model_list list of Model objects. Assess these models.
        #' @param file string or NULL. If not `NULL`, store the resulting plot in this file.
        #' @param fellow_csv logical. Whether to also store the plotted data in a csv.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @export
        assess_center = function(
            data,
            model_list,
            file = NULL,
            fellow_csv = FALSE,
            quiet = FALSE
        )
            ass2d_assess_center(self, private, data, model_list, file, fellow_csv, quiet)
    ),

    private = list(
        # Part of assess method
        get_2d_metric = function(
            data,
            model,
            quiet = FALSE
        )
            ass2d_get_2d_metric(self, private, data, model, quiet)
    )      
)
