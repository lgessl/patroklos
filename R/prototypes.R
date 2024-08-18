#' @title Function interface of a fitter
#' @description This function is the prototype of a fitter with the minimal requirements in its 
#' parameters and return value to work as the `fitter` attribute of a `Model` object.
#' @param x Named numeric matrix. Predictor matrix without `NA`s. Samples correspond to rows.
#' Discrete features are encoded as binary dummy variables.
#' @param y Named list with the response in thee formats:
#' * `"bin"`, a named numeric one-column matrix, binary response,
#' * `"cox"`, a named numeric two-column matrix, with columns `"time_to_event"` and `"event"` 
#' (0 = censoring, 1 = event), the response in the Cox format,
#' * `"true"`, a named numeric one-column matrix, true binary response.
#' 
#' The rownames of `y[["bin"]]` and `y[["true"]]` are a subset of the rownames of `x` and, in 
#' general, do not coincide. Use `intersect_by_names()` to get equal rownames.
#' @param val_error_fun Function to calculate the error of validated predictions. For its interface, 
#' see [`val_error_fun_prototype()`]
#' @param ... Further, fitter-specific hyperparameters.
#' @return An S3 object with underlying class `list`, which we call `fit_obj`. The named list must 
#' have the following element:
#' * `"val_predict"`: a numeric one-column matrix with row names, the validated predictions of the 
#'  (picked) fitted model. The row names are a subset of the row names of `x`. The fitter may 
#' tune hyperparameters and therefore fit multiple models. The `val_predict` attribute must contain
#' of the best validated model among them.
#' 
#' There must be a `predict()` method for the `fit_obj`. See [`predict_generic_prototype()`] for 
#' its interface.
#' @export
fitter_prototype <- function(x, y, val_error_fun, ...){
    stop("This is a prototype function. Do not call it.")
}

#' @title Function interface of a validation error function
#' @description This function is the prototype of a function to calculate the error of validated
#' predictions. It fulfills the minimal requirements in its parameters and return value to work as
#' the `val_error_fun` attribute of a `Model` object.
#' @param y numeric vector. True outcomes.
#' @param y_hat numeric vector. Predicted outcomes. Must have the same length as `y`. Samples in 
#' `y` and `y_hat` must coincide by order.
#' @return A numeric scalar.
val_error_fun_prototype <- function(y, y_hat){
    stop("This is a prototype function. Do not call it.")
}

#' @title Function interface of the S3 generic `predict()` for a `fit_obj`
#' @description This function is the prototype of a `predict()` S3 method for the S3 object 
#' `fit_obj` of a `Model` and lays out the minimal requirements for its parameters and return value.
#' @param object The fitted S3 object.
#' @param newx Numeric matrix. The predictor matrix without `NA`s. Samples correspond to rows. 
#' Discrete features are encoded as binary dummy variables. Column names are the same as in the
#' predictor matrix used for fitting.
#' @param ... Further parameters specific for the S3 class of `object`. Ignored by 
#' `Model$predict()`.
#' @return A numeric one-column matrix with as many rows as `newx`, the predictions of the (picked) 
#' fitted model for the samples in `newx`. The fitter that outputs `object` may fit multiple models
#' and store them in `object`, but the `predict()` method must return the predictions of the best
#' model among them.
predict_method_prototype <- function(object, newx, ...){
    stop("This is a prototype function. Do not call it.")
}

#' @title Function interface for the return value of [`multitune()`]
#' @inheritParams fitter_prototype
#' @param ... atomic vectors. Every vector corresponds to a hyperparameter and holds candidate 
#' values for it. For every combination of hyperparameters, fit a model by calling the `fitter` 
#' parameter of `multitune()`.
#' @return A `multitune_obj` S3 object with underlying class `list`. It fulfills the requirements 
#' of `fitter_prototype()` and, in addition, it has the following elements:
#' * `val_predict_list`: a list of row-named one-column matrices, the validated predictions of the 
#' fitted models.
#' * `lambda`: a character vector, the hyperparameter combinations as a string.
#' * `lambda_min_index`: an integer, the index of the hyperparameter combination of the model with 
#' the lowest validation error.
#' * `lambda_min`: a character, the hyperparameter combination of the model with the lowest
#' validation error.
#' * `val_error`: a numeric vector, the validation errors of the fitted models.
#' * `min_error`: a numeric scalar, the validation error of the model with the lowest validation
#' error.
#' * `fit_obj_list`: a list of `fit_obj`s of the fitted models. If `select = TRUE` in 
#' [`multitune()`], all but the model with the lowest validation error are `NA`.
multitune_output_prototype <- function(x, y, val_error_fun, ...){
    stop("This is a prototype function. Do not call it.")
}