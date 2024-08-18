#' @title An R6 class to assess a model with scalar metrics
#' @description Assess how well a model can predict time to event less than a certain 
#' threshold with a *scalar* metric.
#' @export
AssScalar <- R6::R6Class("AssScalar",
    public = list(

        #' @field metrics Assess the model for these scalar metrics. Check out the 
        #' initializer for possible choices.
        metrics = NULL,
        #' @field prev_range For metrics that need thresholding only consider 
        #' thresholds that yield a prevalence in this range.
        prev_range = NULL,
        #' @field confidence_level Confidence level gamma, e.g. for confidence intervals.
        confidence_level = NULL,
        #' @field benchmark Name and pivot time cutoff of the benchmark `Model`.
        benchmark = NULL,
        #' @field round_digits Round the results in tables to round_digits digits 
        #' after the point.
        round_digits = NULL,
        #' @field file Save the resulting tibble to this csv file. 
        file = NULL,

        #' @description Construct an `AssScalar` R6 object.
        #' @param metrics character. Assess the model for these metrics. For 
        #' currently offered choices see "Usage". If you have a model with 
        #' non-binary output (like the linear predictor of a Cox model), we choose a threshold 
        #' by maximizing the left-most metric in `metrics` *that is made for 
        #' classifiers with binary output* (e.g. precision within `prev_range` below). 
        #' If this cannot be done reasonably, throw an error. Make sure that `hr` precedes 
        #' `hr_ci_ll`, `hr_ci_ul` and `hr_p` in `metrics`; `precision_ci_ll` must precede 
        #' `precision_ci_ul`.
        #' @param prev_range numeric numeric vector of length 2. For metrics that
        #' need thresholding only consider thresholds that yield a prevalence in
        #' this range.
        #' @param confidence_level numeric. The confidence level gamma (e.g. for confidence 
        #' intervals).
        #' @param benchmark list or `NULL`. If not NULL, it is a list with names 
        #' 
        #' * `"name"`: the `name` attribute of the benchmark `Model` in the `model_list` parameter 
        #' of the `assess()` and `assess_center()` method,
        #' * `"prev_range"`: An extra value for the `prev_range` attribute used for the 
        #' benchmark `Model`. Often, we need a higher prevalence for our, new models to gain 
        #' statistical power and be able to significantly outperform the benchmark.
        #' @param file string or NULL. If not `NULL`, save the resulting tibble to this csv file.
        #' @param round_digits numeric. The number of digits to round the results to. 
        #' @return A new `AssScalar` object.
        initialize = function(
            metrics = c("auc", "accuracy", "precision", "prevalence", "precision_ci_ll", 
                "precision_ci_ul", "hr", "hr_ci_ll", "hr_ci_ul", "hr_p", "n_true", "perc_true", 
                "n_samples", "logrank", "threshold"),
            prev_range = c(0, 1),
            confidence_level = 0.95,
            benchmark = NULL,
            file = NULL,
            round_digits = 3
        )
            ass_scalar_initialize(self, private, metrics, prev_range, confidence_level, benchmark, 
                file, round_digits), 

        #' @description Assess a *single* model. 
        #' @param data Data object. Assess on this data. Data must already be read in 
        #' and its `cohort` attribute set.
        #' @param model Model object. Assess this model.
        #' @param quiet logical. Whether to suppress messages.
        #' @return named numeric vector. The calculated metrics.
        assess = function(
            data,
            model,
            quiet = FALSE
        )
            ass_scalar_assess(self, private, data, model, quiet),

        #' @description Wrap `assess()` to assess *multiple* models and store the result.
        #' @param data Data object. Assess on this data. The `cohort` attribute of `data` must be 
        #' set.
        #' @param model_list list of Model objects. Assess these models.
        #' @param quiet logical. Whether to suppress messages.
        #' @return A tibble of shape `(length(model_list) x length(metrics))`.
        #' @export
        assess_center = function(
            data,
            model_list,
            quiet = FALSE
        )
            ass_scalar_assess_center(self, private, data, model_list, 
                mirror, quiet)
    )
)