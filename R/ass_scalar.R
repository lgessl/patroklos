#' @title An R6 class to assess a model with a scalar metric
#' @description Assess how well a model can predict time to event less than a certain 
#' threshold with a *scalar* metric (like the ROC-AUC).
#' @export
AssScalar <- R6::R6Class("AssScalar",
    public = list(

        #' @field metrics Assess the model for these scalar metrics. Check out the 
        #' initializer for possible choices.
        metrics = NULL,
        #' @field prev_range For metrics that need thresholding only consider 
        #' thresholds that yield a prevalence in this range.
        prev_range = NULL,
        #' @field confidence_level Confidence level alpha.
        confidence_level = NULL,
        #' @field benchmark Incorporate the benchmark into the assessment.
        benchmark = NULL,
        #' @field round_digits Round the results in tables to round_digits digits 
        #' after the point.
        round_digits = NULL,
        #' @field file Save the resulting tibble to this csv file. 
        file = NULL,

        #' @description Create a new AssScalar instance.
        #' @param metrics character. Assess the model for these metrics. For 
        #' currently offered choices see "Usage" above. If you have a model with 
        #' non-binary output (like a logistic regression), we choose a threshold 
        #' by maximizing the left-most metric in `metrics` *that is made for 
        #' classifiers with binary output* (e.g. precision within prev_range below). 
        #' If this cannot be done reasonably, throw an error. Make sure that `hr` precedes 
        #' `hr_ci_ll`, `hr_ci_ul` and `hr_p` in `metrics`.
        #' @param prev_range numeric numeric vector of length 2. For metrics that
        #' need thresholding only consider thresholds that yield a prevalence in
        #' this range.
        #' @param confidence_level numeric. The confidence level alpha (e.g. for confidence 
        #' intervals).
        #' @param benchmark character or NULL. If not NULL, include `benchmark` 
        #' (the name of column in the pheno data) in the assessment. Not 
        #' implemented yet.
        #' @param file string or NULL. The name of the csv file to save the 
        #' results to for the *train* cohort.
        #' the results are not saved.
        #' @param round_digits numeric. The number of digits to round the results to. Default is `3`.
        #' @return A new AssScalar object.
        initialize = function(
            metrics = c("auc", "accuracy", "precision", "prevalence", "precision_ci_ll", 
                "HR | lower CI | upper CI | p", "n_true", "perc_true", "n_samples", "logrank", 
                "threshold"),
            prev_range = c(0, 1),
            confidence_level = 0.95,
            benchmark = NULL,
            file = NULL,
            round_digits = 3
        )
            ass_scalar_initialize(self, private, metrics, prev_range, confidence_level, benchmark, 
                file, round_digits), 

        #' @description Assess a *single* model on a data set.
        #' @param data Data object. Assess on this data. Data must already be read in 
        #' and `cohort` attribute set.
        #' @param model Model object. Assess this model.
        #' @param quiet logical. Whether to suppress messages.
        #' @return numeric matrix. The rows correspond to the `model` and optionally 
        #' the benchmark, the columns correspond to `self$metrics`.
        #' @details The AssScalar S3 class is tailored for this function.
        assess = function(
            data,
            model,
            quiet = FALSE
        )
            ass_scalar_assess(self, private, data, model, quiet),

        #' @description Assess *multiple* models on a data set.
        #' @param data Data object. Assess on this data.
        #' If `data$cohort` is `NULL`, assess on the test cohort.
        #' @param model_list list of Model objects. Assess these models.
        #' @param mirror character vector of length 2. If you want to store the 
        #' resulting tibbles, get the test-cohort tibble's file name by mirroring 
        #' `ass_scalar$file` according to `model_tree_mirror`, i.e. replacing 
        #' `model_tree_mirror[1]` by `model_tree_mirror[2]` in the `file`
        #' attribute.
        #' @param quiet logical. Whether to suppress messages.
        #' @return A tibble. Every row stands for one model. 
        #' * If `length(self$metrics) == 1`, the columns hold statistics on the 
        #' single metric (like mean, standard deviation).
        #' * Otherwise columns correspond to `self$metrics`. 
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