#' @title An R6 class for a model
#' @description A Model specifies how a model looks like, prepares the data, fits it, 
#' stores it and predicts from it.
Model <- R6::R6Class("Model",
    public = list(

        #' @field name A telling name for the model.
        name = NULL,
        #' @field directory Store/find the models in this directory.
        directory = NULL,
        #' @field fitter The model fitting function to be used.
        fitter = NULL,
        #' @field split_index Split the given data into training and test cohort
        #' `length(split_index)` times.
        split_index = NULL,
        #' @field time_cutoffs Threshold and censor the outcome accordingly.
        time_cutoffs = NULL,
        #' @field val_error_fun Calculate the error of independently validated predictions.
        val_error_fun = NULL,
        #' @field hyperparams Optional arguments passed to `fitter`.
        hyperparams = NULL,
        #' @field response_colnames Use as column names for the response matrix.
        response_colnames = NULL,
        #' @field include_from_continuous_pheno The names of the continuous variables in the
        #' pheno data (to be) included in the predictor matrix.
        include_from_continuous_pheno = NULL,
        #' @field include_from_discrete_pheno The names of the discrete variables in the
        #' pheno data (to be) included in the predictor matrix.
        include_from_discrete_pheno = NULL,
        #' @field include_expr Whether to include the expression data in the predictor matrix. 
        include_expr = NULL,
        #' @field li_var_suffix Append this to the names of features from the pheno data
        #' when adding them to the predictor matrix.
        li_var_suffix = NULL,
        #' @field create_directory Whether to create `directory` if it does not exist, yet.
        create_directory = NULL,
        #' @field fit_file Store this Model object under this name in `directory`.
        fit_file = NULL,
        #' @field fits A list holding fits (something returned by a fitter like 
        #' a `ptk_zerosum` S3 object).
        fits = NULL,
        #' @field continuous_output Whether the output of the model is continous or binary.
        continuous_output = NULL,
        #' @field combine_n_max_categorical_features Maximum number of categorical features 
        #' to combine 
        combine_n_max_categorical_features = NULL,
        #' @field combined_feature_min_positive_ratio Minimum ratio of positive
        #' observations in a combined (categorical) feature
        combined_feature_min_positive_ratio = NULL,

        #' @description Create a new Model instance.
        #' @param name string. A telling name for the model.
        #' @param fitter function. The model fitting function to be used. Must take `x` and
        #' `y` as first two positional arguments. Further arguments can be passed via
        #' `hyperparams` (below). Its return value must be an S3 object with a `plot()` 
        #' method, and (ideally, for assessment) with a `predict()` method. Default is `NULL`.
        #' @param directory string. The directory to store the models in. For every value in 
        #' `time_cutoffs`, find the corresponding model in a subdirectory named after this value. 
        #' @param split_index integer vector. Split the given data into training and test samples 
        #' `length(split_index)` times, i.e., every index in `split_index` will get its own split. 
        #' @param time_cutoffs numeric vector. This governs how the provided 
        #' response looks like. For every value of `time_cutoffs`, specify a model 
        #' on the following response data:
        #' * For Cox response, censor all patients where the event ouccured after 
        #' `time_cutoffs` at this value and train the specified model.
        #' * For binary response, binarize the outcome depending on whether it 
        #' occured before or after this value.
        #' 
        #' When fitting models to the data, you find every of the 
        #' `length(time_cutoffs)` models in a subdirectory of `directory` named 
        #' after the respective value in `time_cutoffs`.
        #' @param val_error_fun Function used to calculate the error of independently 
        #' validated predictions. Must take two numeric vector of equal length:
        #' `y` and `y_hat`, the true and predicted outcomes, respectively, and 
        #' return a numeric scalar; the lower, the better the model. See
        #' [`error_rate()`] or [`neg_roc_auc()`] for examples. 
        #' @param hyperparams list. Optional arguments passed to `fitter`, e.g. alpha 
        #' in case of an elastic net. Default is `list()`, i.e., no arguments other than `x`, `y`
        #' passed to `fitter`.
        #' @param response_colnames string vector of length 2. Column names for 
        #' the Cox response matrix.
        #' * The first element is the name of the column holding the time until the event or 
        #' censoring, and 
        #' * the second one is the anme of the column holding the event status (1 = event, 0 =
        #' censoring). 
        #' @param include_from_continuous_pheno vector of strings. The names of the 
        #' *continuous* variables in the pheno data (to be) included in the predictor matrix. The
        #' values will be coerced to numeric. Default is `NULL`, which means no continuous pheno
        #' variables are or will be included.
        #' @param include_from_discrete_pheno vector of strings. The names of the *discrete*
        #' variables in the pheno data (to be) included in the predictor matrix. A discrete
        #' variable with n levels will be converted to n-1 binary variables. Default is `NULL`,
        #' which means no discrete pheno variables are or will be included.
        #' @param include_expr logical. Whether to include the expression data in the predictor
        #' matrix.
        #' @param li_var_suffix string. Append this to the names of features from the pheno
        #' data when adding them to the predictor matrix. Default is `"++"`.
        #' @param create_directory logical. Whether to create `directory` if it does not exist, yet. 
        #' Default is `TRUE`.
        #' @param fit_file string. The name of the model-fits file inside `directory`.
        #' Default is `"fit_obj.rds"`.
        #' @param continuous_output logical or NULL. Whether the output of the model 
        #' is continous (like the conditional probabilities of a logistic regression) or
        #' binary (like with a random forest). This piece of information is helpful 
        #' when it comes to assessing a model (Ass2d$assess_center() uses it, e.g.).
        #' @param combine_n_max_categorical_features integer. Maximum number of categorical features
        #' to combine in predicting features.
        #' @param combined_feature_min_positive_ratio numeric. Minimum ratio of positive
        #' observations in a combined (categorical) feature. This attribute together 
        #' with `combine_n_max_categorical_features` governs which combinatorial 
        #' features the predicitor matrix will contain.
        #' @return A `Model` R6 object.
        #' @details Strictly speaking, one `Model` instance specifies
        #' `length(time_cutoffs) * length(split_index)` models. In terms of storing 
        #' and assessing models, we consider the models obtained via repeated 
        #' splitting according to `split_index` as one model; we view models 
        #' obtained via different values of `time_cutoffs`, in contrast, as 
        #' different models; e.g., we can compare them against one another in an 
        #' assessment.
        initialize = function(
            name,
            fitter,
            directory,
            split_index,
            time_cutoffs,
            val_error_fun,
            hyperparams = NULL,
            response_colnames = c("time", "status"),
            include_from_continuous_pheno = NULL,
            include_from_discrete_pheno = NULL,
            include_expr = TRUE,
            li_var_suffix = "++",
            create_directory = TRUE,
            fit_file = "models.rds",
            continuous_output = NULL,
            combine_n_max_categorical_features = 1L,
            combined_feature_min_positive_ratio = 0.04
        )
            model_initialize(self, private, name, fitter, directory, split_index, 
                time_cutoffs, val_error_fun, hyperparams, response_colnames, 
                include_from_continuous_pheno, include_from_discrete_pheno, 
                include_expr, li_var_suffix, create_directory, 
                fit_file, continuous_output,
                combine_n_max_categorical_features, 
                combined_feature_min_positive_ratio),  

        #' @description Fit the model to a data set for all splits into training 
        #' and test cohort. 
        #' @param data Data object. Read it in if needed.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @param msg_prefix string. Prefix for messages. Default is `""`.
        #' @return The `Model` object itself with the `fits` attribute set to a 
        #' list holding the object returned by the `fitter` attribute for every
        #' split. 
        #' @seealso [`training_camp()`].
        fit = function(
            data,
            quiet = FALSE,
            msg_prefix = ""
        )
            model_fit(self, private, data, quiet, msg_prefix),

        #' @description Predict for a data set and all splits into training and 
        #' test cohort. We don't support multiple time cutoffs here. Additonally 
        #' return the true values of the response and, if the `benchmark_col` 
        #' attribute of the `Data` object is not `NULL`, the values of the 
        #' benchmark.
        #' @param data Data object. Specifications on the data. Read it in if 
        #' needed.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @return A list holding:
        #' 
        #' * `"predicted"`: a list of named numeric vectors, the scores output by the model for 
        #'  each split (split index corresponding to list index).
        #' *  "actual": a list of named numeric vectors, for each split the actual values of 
        #' whether time to event was above or below `model$time_cutoffs`, encoded as 1 
        #'  ("high risk") and 0 ("low risk"), respectively. 
        #' * "benchmark": A list of named numeric vectors, for each split the values of the
        #'  benchmark classifier. If `data$benchmark` is NULL, it is an empty list.
        #' 
        #' For every split, the names of all three vectors match.
        #' @importFrom stats predict
        #' @importFrom stats coef
        #' @export
        predict = function(
            data,
            quiet = FALSE
        )
            model_predict(self, private, data, quiet)
    ),

    private = list(
        dummy = "dummy"
    )
)