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
    append_to_includes
){
    # Split into early and late part
    early_bool <- vapply(
        colnames(x),
        function(s) 
            stringr::str_sub(s, -nchar(append_to_includes)) != 
            append_to_includes,
        logical(1)
    ) 
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
    fits <- vector("list", length(fit$lambda) * nrow(hyperparams2))
    dim(fits) <- c(length(fit$lambda), nrow(hyperparams2))
    best_acc <- 0
    best_idx <- c(1L, 1L)
    best_hyperparams <- NULL 
    for (l in seq_along(fit$lambda)) {
        x_late <- cbind(fit$cv_predict[[l]], x[, !early_bool])
        for (m in seq_len(nrow(hyperparams2))) {
           fits[[l, m]] <- do.call(
                fitter2,
                c(list(x = x_late, y = y), as.list(hyperparams2[m, ]))
           )
           acc <- 1-fits[[l, m]]$prediction.error 
           if (acc > best_acc) {
               best_acc <- acc 
               best_idx <- c(l, m)
               best_hyperparams <- c(
                    list("lambda" = fit$lambda[l]), 
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
   stopifnot("predict" %in% sloop::s3_methods_class(class(model1)[1])) 
   stopifnot("predict" %in% sloop::s3_methods_class(class(model2)[1]))
   stopifnot(is.list(best_hyperparams))
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
predict.nested_fit <- function(
    object = NULL,
    newx = NULL,
    ...
){
    stopifnot(inherits(newx, "matrix"))
    stopifnot(!is.null(colnames(newx)))
    stopifnot(is.character(append_to_includes))
    early_bool <- !stringr::str_detect(
        colnames(newx), 
        paste0(append_to_includes, "$")
    )
    x_early <- newx[, early_bool]
    x_late <- cbind(
        predict(object$model1, newx = x_early, s = object$best_hyperparams[1]),
        newx[, !early_bool]
    )
    y <- predict(object$model2, newx = x_late)
    rownames(y) <- rownames(newx)
}