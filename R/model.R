#' @title An R6 class for a model
#' @description A Model specifies how a model looks like, fits and validates it, tunes 
#' hyperparameters, stores it and predicts from it.
Model <- R6::R6Class("Model",
    public = list(

        #' @field name A telling name for the model.
        name = NULL,
        #' @field directory Store/find the `Model` with `fit_obj` set in this directory.
        directory = NULL,
        #' @field fitter Fit and validate the model with this fitting function.
        fitter = NULL,
        #' @field time_cutoffs Threshold and censor the outcome accordingly.
        time_cutoffs = NULL,
        #' @field val_error_fun Calculates the error of the validated predictions.
        val_error_fun = NULL,
        #' @field hyperparams Optional arguments passed to `fitter`.
        hyperparams = NULL,
        #' @field include_from_continuous_pheno The names of the continuous variables in the
        #' pheno data (to be) included in the predictor matrix.
        include_from_continuous_pheno = NULL,
        #' @field include_from_discrete_pheno The names of the discrete variables in the
        #' pheno data (to be) included in the predictor matrix.
        include_from_discrete_pheno = NULL,
        #' @field include_expr Whether to include the expression data in the predictor matrix. 
        include_expr = NULL,
        #' @field combine_n_max_categorical_features Maximum number of categorical features 
        #' to combine.
        combine_n_max_categorical_features = NULL,
        #' @field combined_feature_min_positive_ratio Minimum ratio of positive
        #' observations in a combined (categorical) feature.
        combined_feature_min_positive_ratio = NULL,
        #' @field enable_imputation Overrides the `imputer` attribute of the `Data` object.
        enable_imputation = NULL,
        #' @field fit_obj The fitted object, something returned by a fitter like a `ptk_zerosum` S3 
        #' object.
        fit_obj = NULL,
        #' @field create_directory Whether to create `directory` if it does not exist, yet.
        create_directory = NULL,
        #' @field file Store this Model object under this name in `directory`.
        file = NULL,
        #' @field li_var_suffix Append this to the names of features from the pheno data
        #' when adding them to the predictor matrix.
        li_var_suffix = NULL,

        #' @description Create a new Model instance.
        #' @param name string. A telling name for the model.
        #' @param fitter function. Fit the model to the data and validate it with this fitting 
        #' function. See [`fitter_prototype()`] for its required interface. patroklos provides 
        #' two fitters out of the box: [`ptk_zerosum()`], a wrapper around [`zeroSum::zeroSum`] and 
        #' [`ptk_ranger()`], a wrapper around [`ranger::ranger()`]. To tune more than just one 
        #' combination of hyperparameters, decorate a fitter with [`multitune()`].
        #' @param directory string. Store/find the `Model` with `fit_obj` set in this directory.
        #' @param time_cutoffs numeric vector. A model-agnostic hyperparameter that changes the 
        #' response during *training* as follows: For every value of `time_cutoffs`, specify a 
        #' model on the following response data:
        #' * For Cox response, censor all patients where the event occurred after 
        #' `time_cutoffs` at this value and train the specified model.
        #' * For binary response, binarize the outcome depending on whether it 
        #' occurred before or after this value.
        #' 
        #' We already tune this hyperparameter and only store the best model according to 
        #' validation as `fit_obj` and report the chosen time cutoff as `time_cutoff` attribute 
        #' in it.
        #' @param val_error_fun Function to calculate the error of validated predictions. For its 
        #' interface, see [`val_error_fun_prototype()`].
        #' @param hyperparams list. Optional arguments passed to `fitter`, e.g. alpha 
        #' in case of an elastic net.
        #' @param include_from_continuous_pheno vector of strings. The names of the 
        #' *continuous* variables in the pheno data (to be) included in the predictor matrix. 
        #' Default is `NULL`, which means no continuous pheno variables are or will be included.
        #' @param include_from_discrete_pheno vector of strings. The names of the *discrete*
        #' variables in the pheno data (to be) included in the predictor matrix. A discrete
        #' variable with n levels will be dichotomized into n-1 binary dummy variables. Default is 
        #' `NULL`, which means no discrete pheno variables are or will be included.
        #' @param include_expr logical. Whether to include the expression data in the predictor
        #' matrix.
        #' @param combine_n_max_categorical_features integer. Maximum number of categorical features
        #' to combine in predicting features.
        #' @param combined_feature_min_positive_ratio numeric. Minimum ratio of positive
        #' observations in a combined (categorical) feature. This attribute together 
        #' with `combine_n_max_categorical_features` governs which combined categorical
        #' features the predictor matrix will contain: add a combination of the levels of distinct 
        #' categorical features to the predictor matrix after imputation if at most 
        #' `combine_n_max_categorical_features` are involved in the combination and if the 
        #' combination is (expected to be) there in at least `combined_feature_min_positive_ratio`
        #' of the samples.
        #' @param enable_imputation logical. If `FALSE`, it overrides the `imputer` attribute of 
        #' the `Data` object and we do not impute.
        #' @param file string. The name of the model-fit_obj file inside `directory`.
        #' Default is `"fit_obj.rds"`.
        #' @param create_directory logical. Whether to create `directory` if it does not exist, yet. 
        #' Default is `TRUE`.
        #' @param li_var_suffix string. Append this to the names of features from the pheno
        #' data when adding them to the predictor matrix. Default is `"++"`.
        #' @return A `Model` R6 object.
        initialize = function(
            name,
            fitter,
            directory,
            time_cutoffs,
            val_error_fun,
            hyperparams = NULL,
            include_from_continuous_pheno = NULL,
            include_from_discrete_pheno = NULL,
            include_expr = TRUE,
            combine_n_max_categorical_features = 1L,
            combined_feature_min_positive_ratio = 0.04,
            enable_imputation = TRUE,
            file = "model.rds",
            create_directory = TRUE,
            li_var_suffix = "++"
        )
            model_initialize(self, private, name, fitter, directory, time_cutoffs, 
                val_error_fun, hyperparams, include_from_continuous_pheno, 
                include_from_discrete_pheno, include_expr, li_var_suffix, 
                create_directory, file, 
                combine_n_max_categorical_features, 
                combined_feature_min_positive_ratio, enable_imputation),  

        #' @description Fit the model to a data set, validate it and tune hyperparameters.
        #' @param data Data object. Read it in if needed.
        #' @param update_model_shell logical. If `TRUE` and we find a stored model with `fit_obj` 
        #' not being NULL, we set the `fit_obj` attribute of the model to the found `fit_obj`
        #' and save it. This way, we can keep stored `Model`s up-to-date with changes in the 
        #' `Model` class.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @param msg_prefix string. Prefix for messages. Default is `""`.
        #' @return The `Model` object itself with the `fit_obj` attribute set 
        #' to the object tuned over `time_cutoffs`, `combine_n_max_categorical_features` 
        #' and `hyperparams`.
        #' @seealso [`training_camp()`].
        fit = function(
            data,
            update_model_shell = FALSE,
            quiet = FALSE,
            msg_prefix = ""
        )
            model_fit(self = self, private = private, data = data, 
                update_model_shell = update_model_shell, quiet = quiet, msg_prefix = msg_prefix),

        #' @description Predict for a data set. 
        #' @param data Data object. Specifications on the data. Read it in if 
        #' needed.
        #' @param quiet logical. Whether to suppress messages. Default is `FALSE`.
        #' @return A named list of length 3: 
        #' * `"predicted"`: named numeric vector with predicted response,
        #' * `"actual"`: named numeric vector with actual response,
        #' * `"cox_mat"`: named matrix with columns `"time_to_event"`, `"event"`, `"hazard"`, 
        #' where hazard again is the predicted response. This matrix is helpful to calculate 
        #' hazard ratios or logrank p-values.
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