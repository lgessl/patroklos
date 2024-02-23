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
        #' @param file string or NULL. The name of the csv file to save the results to. If `NULL`,
        #' the results are not saved. If you supply an AssScalar object to `AssScalar_assess_center()`,
        #' specify the path for the *train* cohort.
        #' @param round_digits numeric. The number of digits to round the results to. Default is `3`.
        #' @return A new AssScalar object.
        initialize = function(
            metric,
            pivot_time_cutoff,
            lambda = "lambda.min",
            benchmark = NULL,
            file = NULL,
            round_digits = 3
        ){
            stopifnot(is.character(metric))
            stopifnot(is.numeric(pivot_time_cutoff))
            stopifnot(pivot_time_cutoff > 0)
            stopifnot(is.character(lambda) || is.numeric(lambda))
            stopifnot(is.null(benchmark) || is.character(benchmark))
            stopifnot(is.null(file) || is.character(file))
            stopifnot(is.numeric(round_digits) && round_digits >= 0)
            self$metric <- metric
            self$pivot_time_cutoff <- pivot_time_cutoff
            self$lambda <- lambda
            self$benchmark <- benchmark
            self$file <- file
            self$round_digits <- round_digits
        }, 

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
        ){
            AssScalar_assess(
                ass_scalar = self,
                data = data,
                model = model,
                quiet = quiet
            )
        },

        #' @description Assess *multiple* models (with multiple splits) on a data set.
        #' @param ass_scalar AssScalar S3 object. See the constructor [`AssScalar()`] for more 
        #' details.
        #' @param data Data object. Assess on this data.
        #' If `data$cohort` is `NULL`, assess on the test cohort.
        #' @param model_list list of Model objects. Assess these models.
        #' @param model_tree_mirror character vector of length 2. If you want to store the 
        #' resulting tibbles, get the test-cohort tibble's file name by mirroring 
        #' `ass_scalar$file` according to `model_tree_mirror` (see [`mirror_path()`]).
        #' @param quiet logical. Whether to suppress messages.
        #' @return A tibble. Every row holds the metric (and more analysis across the splits 
        #' like standard deviation, minimum, maximum value) for one model in `model_list` 
        #' and every time cutoff specified for this model.
        #' @export
        assess_center = function(
            data,
            model_list,
            model_tree_mirror = c("models", "results"),
            quiet = FALSE
        ){
            AssScalar_assess_center(
                ass_scalar = self,
                data = data,
                model_list = model_list,
                model_tree_mirror = model_tree_mirror,
                quiet = quiet
            )
        }
    )
)