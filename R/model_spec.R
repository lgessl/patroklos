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
        #' @field optional_fitter_args Optional arguments passed to `fitter`.
        optional_fitter_args = NULL,
        #' @field response_type The type of response to be used.
        response_type = NULL,
        #' @field response_colnames Use as column names for the response matrix.
        response_colnames = NULL,
        #' @field include_from_continuous_pheno The names of the continuous variables in the
        #' pheno data (to be) included in the predictor matrix.
        include_from_continuous_pheno = NULL,
        #' @field include_from_discrete_pheno The names of the discrete variables in the
        #' pheno data (to be) included in the predictor matrix.
        include_from_discrete_pheno = NULL,
        #' @field append_to_includes Append this to the names of features from the pheno data
        #' when adding them to the predictor matrix.
        append_to_includes = NULL,
        #' @field create_directory Whether to create `directory` if it does not exist, yet.
        create_directory = NULL,
        #' @field plot_file Store the plots resulting from `plot(fit_obj)` in `directory` under
        #' this name.
        plot_file = NULL,
        #' @field plot_ncols Arrange the above mentioned plots in this number of columns. 
        plot_ncols = NULL,
        #' @field plot_title_line Pass this as the `line` argument to [`graphics::title()`]
        #' when calling `plot(fit_obj)`.
        plot_title_line = NULL,
        #' @field fit_file Store this Model object under this name in `directory`.
        fit_file = NULL,

        #' @description Create a new Model instance.
        #' @param name string. A telling name for the model.
        #' @param directory string. The directory to store the models in. For every value in 
        #' `time_cutoffs`, find the corresponding model in a subdirectory named after this value. 
        #' @param fitter function. The model fitting function to be used. Must take `x` and
        #' `y` as first two positional arguments. Further arguments can be passed via
        #' `optional_fitter_args` (below). Its return value must be an S3 object with a `plot()` 
        #' method, and (ideally, for assessment) with a `predict()` method. Default is `NULL`.
        #' @param time_cutoffs numeric vector.
        #' * If `response_type == "survival_censored"`: For every value in `time_cutoffs`, censor
        #' all patients where the event ouccured after `time_cutoffs` at this value and train the 
        #' specified model.
        #' * If `response_type == "binary"`: For every value in `time_cutoffs`, binarize the 
        #' outcome depending on whether it occured before or after this value, and train the 
        #' specified model.
        #' @param split_index integer vector. Split the given data into training and test samples 
        #' `length(split_index)` times, i.e., every index in `split_index` will get its own split. 
        #' @param optional_fitter_args list. Optional arguments passed to `fitter`, e.g. alpha 
        #' in case of an elastic net. Default is `list()`, i.e., no arguments other than `x`, `y`
        #' passed to `fitter`.
        #' @param response_type string. The type of response to be used. One of `"binary"` or
        #' `"survival_censored"`. Default is `NULL`.
        #' @param response_colnames string vector of length 2. If `response_type == "survival_censored"`,
        #' use as column names for the response matrix. 
        #' * The first element is the name of the column holding the time until the event or 
        #' censoring, and 
        #' * the second one is the anme of the column holding the event status (1 = event, 0 =
        #' censoring). 
        #' Default is `c("time", "status")`.
        #' @param include_from_continuous_pheno vector of strings. The names of the 
        #' *continuous* variables in the pheno data (to be) included in the predictor matrix. The
        #' values will be coerced to numeric. Default is `NULL`, which means no continuous pheno
        #' variables are or will be included.
        #' @param include_from_discrete_pheno vector of strings. The names of the *discrete*
        #' variables in the pheno data (to be) included in the predictor matrix. A discrete
        #' variable with n levels will be converted to n-1 binary variables. Default is `NULL`,
        #' which means no discrete pheno variables are or will be included.
        #' @param append_to_includes string. Append this to the names of features from the pheno
        #' data when adding them to the predictor matrix. Default is `"++"`.
        #' @param create_directory logical. Whether to create `directory` if it does not exist, yet. 
        #' Default is `TRUE`.
        #' @param plot_file string. Store the plot resulting from `plot(fit_obj)` in `directory`
        #' under this name. Default is `"training_error.pdf"`.
        #' @param plot_ncols integer. The number of columns in the plot. Default is `2`.
        #' @param plot_title_line numeric or NULL. Pass this as the `line` argument to [`graphics::title()`] 
        #' after calling `plot(fit_obj)`. This is the distance (in inches) between the title text and 
        #' the upper limit of the figure. Default is `2.5`.
        #' @param fit_file string. The name of the model-fits file inside `directory`.
        #' Default is `"fit_obj.rds"`.
        #' @return A ModelSpec S3 object.
        #' @details Strictly speaking, one `Model` instance specifies
        #' `length(time_cutoffs) * length(split_index)` models. In terms of storing and assessing models,
        #' we consider the models obtained via repeated splitting according to `split_index` as one 
        #' model; repeated splitting serves the purpose of getting more reliable estimates of its
        #' performance. We view models obtained via different values of `time_cutoffs`, in 
        #' contrast, as different models; e.g., we can compare them against one another in an 
        #' assessment.
        initialize = function(
            name,
            fitter,
            directory,
            split_index,
            time_cutoffs,
            optional_fitter_args = NULL,
            response_type = c("binary", "survival_censored"),
            response_colnames = c("time", "status"),
            include_from_continuous_pheno = NULL,
            include_from_discrete_pheno = NULL,
            append_to_includes = "++",
            create_directory = TRUE,
            plot_file = "training_error.pdf",
            plot_ncols = 2,
            plot_title_line = 2.5,
            fit_file = "fit_obj.rds"
        ){
            response_type <- match.arg(response_type)
            stopifnot(all(time_cutoffs >= 0))
            stopifnot(split_index >= 1)
            stopifnot(is.character(name))
            check_fitter(fitter, optional_fitter_args)
            stopifnot(is.numeric(time_cutoffs))
            stopifnot(is.numeric(split_index))
            stopifnot(is.character(response_type))
            stopifnot(is.character(response_colnames))
            stopifnot(is.character(include_from_continuous_pheno) || is.null(include_from_continuous_pheno))
            stopifnot(is.character(include_from_discrete_pheno) || is.null(include_from_discrete_pheno))
            stopifnot(is.character(append_to_includes))
            stopifnot(is.character(directory))
            stopifnot(is.logical(create_directory))
            stopifnot(is.numeric(time_cutoffs) || is.null(time_cutoffs))
            stopifnot(is.character(plot_file))
            stopifnot(is.numeric(plot_ncols))
            stopifnot(is.numeric(plot_title_line) || is.null(plot_title_line))
            stopifnot(is.character(fit_file))

            self$name <- name
            self$directory <- directory
            self$fitter <- fitter
            self$split_index <- split_index
            self$time_cutoffs <- time_cutoffs
            self$optional_fitter_args <- optional_fitter_args
            self$response_type <- response_type
            self$response_colnames <- response_colnames
            self$include_from_continuous_pheno <- include_from_continuous_pheno
            self$include_from_discrete_pheno <- include_from_discrete_pheno
            self$append_to_includes <- append_to_includes
            self$create_directory <- create_directory
            self$plot_file <- plot_file
            self$plot_ncols <- plot_ncols
            self$plot_title_line <- plot_title_line
            self$fit_file <- fit_file
        }
    )
)

at_time_cutoff <- function(
    model,
    time_cutoff
){
    if(!time_cutoff %in% model$time_cutoffs){
        stop("time_cutoff must be one of model$time_cutoffs")
    }
    model_cutoff <- model # Rely on copy-on-modify
    model_cutoff$name <- paste0(model$name, "@", time_cutoff)
    model_cutoff$time_cutoffs <- time_cutoff
    model_cutoff$directory <- file.path(
        model$directory, 
        stringr::str_replace(as.character(time_cutoff), "\\.", "-")
    )
    return(model_cutoff)
}