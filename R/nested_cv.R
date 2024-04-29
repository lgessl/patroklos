#' @title Nested cross-validation for second-stage OOB predictions
#' @description Perform a nested cross-validation for a late-integration scheme,
#' i.e., perform a cross-validation for the early model and then train second-stage
#' models on cross-validated predictions and evaluate the entire models via 
#' out-of-bag (OOB) predictions of the second-stage models (usually a random 
#' forest).
#' @param x A numeric matrix holding the predictor features: rows are samples and
#' columns are features.
#' @param y A numeric vector holding the response variable. 
#' @param fitter1 A function that fits the early model. It must be able to perform 
#' a cross validation and return an S3 class with
#' * a `cv_predict` attribute, a list of vectors, for every lambda the
#' cross-validated predictions.
#' * a `lambda` attribute, a numeric vector holding the lambda values, and 
#' * a `predict` method.
#' @param fitter2 A function that fits the late model. It must return an S3 object
#' with a `predict` method.
#' @param hyperparams1 A named list with hyperpamaters we will pass to `fitter1`
#' to obtain cross-validated predictions from the early model.
#' @param hyperparams2 A named list with hyperparameters for the late model. 
#' Unlike `hyperparams1`, we call `fitter2` for every combination of values in 
#' `hyperparams2` and lambda value from `fitter1`.
#' @param n_folds An integer specifying the number of folds for the cross-validation.
#' @param append_to_includes A character string to filter columns in `x` we 
#' provide to the late model as predictors: all columns with names matching
#' `append_to_includes` as a regular expression.
#' @param pseudo_cv A logical indicating whether we use pseudo cross-validation.
#' I.e., `fitter1` performs one `n_folds` cross validation, get a new training 
#' set by replacing `fitter1`'s input features by the cross-validated predictions 
#' and procede with model fitting and getting out-of-bag predictions. This a is 
#' not real independent testing since the OOB prediction for a given sample 
#' according to `fitter2` involves a model that saw the sample during training.
#' Yet, this pseudo method saves us a factor of `n_folds` in computational time 
#' compared to a classical nested cross-validation. The latter will follow in 
#' future versions soon. 
#' @return An S3 object with class `nested_fit`, the model with the best 
#' performance according to the out-of-bag (OOB) predictions based on cross-validated
#' predictions from the early model.   
#' @export
#' @details This function does hyperparameter tuning for a nested model, i.e.,
#' a so-called early model makes predictions from the high-dimensional part of 
#' data (e.g. RNA-seq, Nanostring). We then provide these predictions as a 
#' one-dimensional continuous feature together with new features to a late model 
#' *that enables a performance estimate via out-of-bag (OOB) predictions 
#' (typically a random forest)*.       
nested_cv_oob <- function(
    x,
    y,
    fitter1,
    fitter2,
    hyperparams1,
    hyperparams2,
    n_folds,
    append_to_includes,
    pseudo_cv = TRUE
){
    # Input checks
    stopifnot(is.matrix(x) && is.numeric(x))
    stopifnot(is.numeric(y) || is.factor(y))
    if (is.vector(y)) {
        stopifnot(length(y) == nrow(x))
    } else if (is.matrix(y)) {
        stopifnot(nrow(y) == nrow(x))
    } else {
        stop("y must be a vector or a matrix.")
    }
    stopifnot(length(unique(y)) != 1)
    stopifnot(is.function(fitter1) && is.function(fitter2))
    stopifnot(is.list(hyperparams1) && is.list(hyperparams2))
    n_folds <- as.integer(n_folds)
    stopifnot(n_folds > 1)
    stopifnot(is.character(append_to_includes))
    stopifnot(is.character(append_to_includes))
    stopifnot(is.logical(pseudo_cv))
    if (!pseudo_cv) {
        stop("Right now, only pseudo cross-validation is supported.")
    }
    early_bool <- get_early_bool(x, append_to_includes) 
    x_early <- x[, early_bool] 
    if(is.null(hyperparams2$expand) || hyperparams2$expand){
        hyperparams2 <- expand.grid(hyperparams2)
    }
    # First stage
    fit <- do.call(
        fitter1,
        c(list(x = x_early, y = y, nFold = n_folds), hyperparams1)
    )
    # Second stage
    n_lambda <- length(fit$lambda)
    n_hyper2 <- nrow(hyperparams2)
    fits <- vector("list", length(fit$lambda) * nrow(hyperparams2))
    best_acc <- 0
    best_idx <- 0
    best_hyperparams <- NULL 
    for (l in seq(n_lambda)) {
        x_late <- cbind(fit$cv_predict[[l]], x[, !early_bool])
        for (m in seq(n_hyper2)) {
            idx <- (l-1)*n_hyper2 + m 
            fits[[idx]] <- do.call(
                fitter2,
                c(list(x = x_late, y = y), as.list(hyperparams2[m, ]))
            ) 
            acc <- mean(fits[[idx]]$predictions == y) 
            if(is.nan(acc))
                stop("The S3 object returned by `fitter2` must have a `predictions` 
                    attribute.")
            if (acc > best_acc) {
                best_acc <- acc 
                best_idx <- idx
                best_hyperparams <- c(
                    list("lambda_index" = l, "lambda" = fit$lambda[l]), 
                    as.list(hyperparams2[m, ])
                )
            }
        } 
    }
    nested_fit(
        model1 = fit, 
        model2 = fits[[best_idx]], 
        best_hyperparams = best_hyperparams, 
        append_to_includes = append_to_includes
    )
}

#' @title Construct a nested_fit S3 object
#' @description Construct a nested_fit object that holds the early and late model
#' plus their hyperparameters.
#' @param model1 An S3 object with a `predict` method. The early model.
#' @param model2 An S3 object with a `predict` method. The late model.
#' @param best_hyperparams A named list with the best hyperparameters for the 
#' late model.
#' @param append_to_includes A character string matching exactly the names of 
#' the features we only feed into the late model.
#' @return An S3 object with class `nested_fit'.
#' @details Nest two models. The first model takes part of the features (those 
#' not matching pattern) to make a prediction that we in turn feed together with 
#' the rest of the features to the second model. This gives us a final prediction.  
#' Both models must be S3 objects with a `predict` method.
#' @export
nested_fit <- function(
    model1,
    model2,
    best_hyperparams,
    append_to_includes
){
    methods1 <- unlist(lapply(class(model1), function(cl) sloop::s3_methods_class(
        cl)$generic)) 
    methods2 <- unlist(lapply(class(model2), function(cl) sloop::s3_methods_class(
        cl)$generic)) 
    if (!("predict" %in% methods1)) {
        stop("`model1` must have a predict method.")
    }
    if (!("predict" %in% methods2)) {
        stop("`model2` must have a predict method.")
    }
    stopifnot(is.list(best_hyperparams))
    stopifnot(is.character(append_to_includes))
    structure(
        list(
            "model1" = model1,
            "model2" = model2,
            "best_hyperparams" = best_hyperparams,
            "append_to_includes" = append_to_includes
        ),
        class = c("nested_fit", "list")
    )
}

#' @title Predict method for nested_fit objects
#' @param object An S3 object with class `nested_fit`. The two-stage nested
#' model.
#' @param newx A numeric matrix with the same number of columns and column names 
#' as the matrix used to fit the model. Rows correspond to samples, columns to 
#' features.
#' @param ... Other arguments for the normal predict function if the fit is not 
#' a nested_fit object.
#' @return A numeric vector with the prediction for every sample.
#' @importFrom stats predict
#' @export  
predict.nested_fit <- function(
    object = NULL,
    newx = NULL,
    ...
){
    stopifnot(inherits(newx, "matrix"))
    stopifnot(!is.null(colnames(newx)))
    early_bool <- get_early_bool(newx, object$append_to_includes)
    x_early <- newx[, early_bool]
    x_late <- cbind(
        predict(object = object$model1, newx = x_early, 
            s = object$best_hyperparams$lambda_index),
        newx[, !early_bool]
    )
    y <- predict(object = object$model2, data = x_late, newx = x_late)
    if (!is.numeric(y)) 
        y <- y[[1]]
    if (!is.numeric(y) || length(y) != nrow(newx))
        stop("The predict method of `model2` must return a numeric vector of same 
            length as the number of rows in `newx` or a list with the first element 
            being the former.")
    dim(y) <- NULL
    names(y) <- rownames(newx)
    y
}
