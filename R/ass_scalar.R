#' @title An R6 class to assess a model with a scalar metric
#' @description Assess how well a model can predict time to event less than a certain 
#' threshold with a *scalar* metric (like the ROC-AUC).
#' @export
AssScalar <- R6::R6Class("AssScalar",
    public = list(

        #' @field metric Scalar metric.
        metric = NULL,
        #' @field pivot_time_cutoff Assess classifying time to event less than 
        #' pivot_time_cutoff.
        pivot_time_cutoff = NULL,
        #' @field lambda Assess the model with regularization parameter lambda.
        lambda = NULL,
        #' @field benchmark Incorporate the benchmark into the assessment.
        benchmark = NULL,
        #' @field round_digits Round the results in tables to round_digits digits 
        #' after the point.
        round_digits = NULL,
        #' @field file Save the resulting tibble to this csv file. This is the 
        #' file for the *train* cohort, we infer the the file name of the test 
        #' cohort by mirroring this file name. See the `assess_center()` method 
        #' for more details.
        file = NULL,

        #' @description Create a new AssScalar instance.
        #' @param metric character. Name of the function used to calculate the metric. It must 
        #' take two numeric vectors of the same length as arguments: 
        #' 1. `predicted`, the predictions by a model, higher values indicate a more positive 
        #' prediction,
        #' 2. `actual`, the true values with the positive class coded as 1, the negative one as 0. 
        #' @param pivot_time_cutoff numeric. Assess for classifying time to event < 
        #' `pivot_time_cutoff`.
        #' @param lambda Assess the model belonging to this `lambda` in the cross validation 
        #' (if there was one).
        #' @param benchmark character or NULL. If not NULL, include `benchmark` (the name of column 
        #' in the pheno data) in the assessment.
        #' @param file string or NULL. The name of the csv file to save the 
        #' results to for the *train* cohort.
        #' the results are not saved.
        #' @param round_digits numeric. The number of digits to round the results to. Default is `3`.
        #' @return A new AssScalar object.
        initialize = function(
            metric,
            pivot_time_cutoff,
            lambda = "lambda.min",
            benchmark = NULL,
            file = NULL,
            round_digits = 3
        )
            ass_scalar_initialize(self, private, metric, pivot_time_cutoff, 
                lambda, benchmark, file, round_digits), 

        #' @description Assess a *single* model (with multiple splits) on a data set.
        #' @param data Data object. Assess on this data. Data must already be read in 
        #' and `cohort` attribute set.
        #' @param model Model object. Assess this model, multiple splits are supported.
        #' @param quiet logical. Whether to suppress messages.
        #' @return numeric vector. For every split index the desired metric.
        #' @details The AssScalar S3 class is tailored for this function.
        assess = function(
            data,
            model,
            quiet = FALSE
        )
            ass_scalar_assess(self, private, data, model, quiet),

        #' @description Assess *multiple* models (with multiple splits) on a data set.
        #' @param data Data object. Assess on this data.
        #' If `data$cohort` is `NULL`, assess on the test cohort.
        #' @param model_list list of Model objects. Assess these models.
        #' @param mirror character vector of length 2. If you want to store the 
        #' resulting tibbles, get the test-cohort tibble's file name by mirroring 
        #' `ass_scalar$file` according to `model_tree_mirror`, i.e. replacing 
        #' `model_tree_mirror[1]` by `model_tree_mirror[2]` in the `file`
        #' attribute.
        #' @param quiet logical. Whether to suppress messages.
        #' @return A tibble. Every row holds the metric (and more analysis across the splits 
        #' like standard deviation, minimum, maximum value) for one model in `model_list` 
        #' and every time cutoff specified for this model.
        #' @export
        assess_center = function(
            data,
            model_list,
            mirror = c("models", "results"),
            quiet = FALSE
        )
            ass_scalar_assess_center(self, private, data, model_list, 
                mirror, quiet)
    )
)