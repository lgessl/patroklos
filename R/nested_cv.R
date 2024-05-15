#' @title Nested cross-validation for second-stage OOB predictions
#' @description Perform a nested cross-validation for a late-integration scheme,
#' i.e., perform a cross-validation for the early model and then train second-stage
#' models on cross-validated predictions and evaluate the entire models via 
#' out-of-bag (OOB) predictions of the second-stage models (usually a random 
#' forest).
#' @param x A numeric matrix holding the predictor features: rows are samples and
#' columns are features.
#' @param y A numeric vector holding the response variable. 
#' @param fitter1 A *patroklos-compliant fitter with integrated CV* (for what this 
#' means, see [`ptk_zerosum()`]) to fit the early model.
#' @param fitter2 A *patroklos-compliant fitter with validated predictions* (
#' for what this means, see [`ptk_ranger()`]) to fit the late model.
#' @param hyperparams1 A named list with hyperparaters we will pass to `fitter1`.
#' @param hyperparams2 A named list with hyperparameters for the late model. 
#' Unlike `hyperparams1`, we call `fitter2` for every combination of values in 
#' `hyperparams2` and lambda value from `fitter1`.
#' @param oob A logical vector of length 2. If the first element is `TRUE`, train 
#' the late model on out-of-bag, else cross-validated predictions from the early
#' model. If the second element is `TRUE`, evaluate the entire model on OOB,
#' else cross-validated predictions.
#' @return An S3 object with class `nested_fit`, the model with the best 
#' performance according to the out-of-bag (OOB) predictions based on cross-validated
#' predictions from the early model.   
#' @export
#' @details This function does hyperparameter tuning for a nested model, i.e.,
#' a so-called early model makes predictions from the high-dimensional part of 
#' data (e.g. RNA-seq, Nanostring), we then provide these predictions as a 
#' one-dimensional feature together with new features to a late model. Both the 
#' early and late model try to predict `y`. To not provide the late model overly 
#' optimistic (since overfitted) predictions during training, we feed its 
#' training algorithm with values comparable to those we would observe for 
#' independent test samples, i.e. either cross-validated or out-of-bage (OOB) 
#' predictions. To evaluate the overall model, we do a second cross-validation 
#' or use OOB predictions. 
#' 
#' Note that the predictions we get for the nested model are not predictions as 
#' one would observe for independent test samples:
#' Let's fix sample i. We get the OOB/CV prediction for sample i from models/
#' a model whose training algorithm didn't see sample i itself. But it probably 
#' saw the prediction for a sample j != i according to an early model whose 
#' training algorithm had seen sample i. Hence the term "pseudo" in the name of 
#' this function. This heuristic saves a factor `n_folds` computation time
#' compared to a full nested cross-validation.  
nested_pseudo_cv <- function(
    x,
    y,
    fitter1,
    fitter2,
    hyperparams1,
    hyperparams2,
    oob = c(FALSE, TRUE)
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
    if (length(unique(y)) != 2) {
        stop("y must be binary. but your y has the unique elements: ", 
            paste(unique(y), collapse = ", "), ". (length(y) = ", length(y), 
            ")")
    }
    stopifnot(is.function(fitter1) && is.function(fitter2))
    stopifnot(is.list(hyperparams1) && is.list(hyperparams2))
    stopifnot(is.logical(oob) && length(oob) == 2)

    val_predict_name <- ifelse(oob, "oob_predict", "cv_predict")
    val_predict_name[1] <- paste0(val_predict_name[1], "_list")
    early_bool <- get_early_bool(x) 
    li_var_suffix <- attr(x, "li_var_suffix")
    x_early <- x[, early_bool]
    attr(x_early, "li_var_suffix") <- li_var_suffix

    # First stage
    fit <- do.call(
        fitter1,
        c(list(x = x_early, y = y), hyperparams1)
    )

    # Second stage
    n_lambda <- length(fit[[val_predict_name[1]]]) # Partition for-loop
    hyperparams <- expand.grid(c(
        hyperparams2, 
        list("lambda_index" = seq(n_lambda))
    ))
    accuracy <- rep(NA, nrow(hyperparams)) 
    fits <- vector("list", nrow(hyperparams))
    n_hyper2 <- as.integer(nrow(hyperparams)/n_lambda)
    n_class_hyper2 <- length(hyperparams2)
    hyperparams[["lambda"]] <- fit$lambda[hyperparams[["lambda_index"]]]
    for (l in seq(n_lambda)) {
        x_late <- cbind(fit[[val_predict_name[1]]][[l]], x[, !early_bool])
        if (ncol(x_late) != sum(!early_bool)+1)
            stop("Something went wrong with adding the early model's predictions.")
        for (m in seq(n_hyper2)) {
            idx <- (l-1)*n_hyper2 + m 
            fits[[idx]] <- do.call(
                fitter2,
                c(
                    list(x = x_late, y = y), 
                    as.list(hyperparams[m, seq(n_class_hyper2)])
                )
            ) 
            acc <- mean(fits[[idx]][[val_predict_name[2]]] == y) 
            if(is.nan(acc))
                stop("The S3 object returned by `fitter2` must have a `predictions` 
                    attribute.")
            accuracy[idx] <- acc 
        }
    } 
    best_idx <- which.max(accuracy)
    cv_oob <- ifelse(oob[2], "oob", "cv")
    hyperparams[[paste0("overall_", cv_oob, "_accuracy")]] <- accuracy

    nested_fit(
        model1 = fit, 
        model2 = fits[[best_idx]], 
        search_grid = hyperparams,
        best_hyperparams = hyperparams[best_idx, ]
    )
}

#' @title Construct a nested_fit S3 object
#' @description Construct a nested_fit object that holds the early and late model
#' plus their hyperparameters.
#' @param model1 An S3 object with a `predict` method. The early model.
#' @param model2 An S3 object with a `predict` method. The late model.
#' @param search_grid A data frame with the hyperparamters for the early
#' (`lambda`) and the late model, i.e., one row corresponds to one trained 
#' model or to one combinations of hyperparameters. Plus an extra column for the
#' cross-validated or out-of-bag (OOB) accuracy.
#' @param best_hyperparams A named list with the best hyperparameters.
#' @return An S3 object with class `nested_fit'.
#' @details Nest two models. The first model takes part of the features (those 
#' not matching pattern) to make a prediction that we in turn feed together with 
#' the rest of the features to the second model. This gives us a final prediction.  
#' Both models must be S3 objects with a `predict` method.
#' @export
nested_fit <- function(
    model1,
    model2,
    search_grid,
    best_hyperparams
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
    stopifnot(is.data.frame(search_grid))
    stopifnot(is.list(best_hyperparams))
    structure(
        list(
            "model1" = model1,
            "model2" = model2,
            "search_grid" = search_grid,
            "best_hyperparams" = best_hyperparams
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
    object,
    newx,
    ...
){
    stopifnot(inherits(newx, "matrix"))
    stopifnot(!is.null(colnames(newx)))
    early_bool <- get_early_bool(newx)
    x_early <- newx[, early_bool]
    attr(x_early, "li_var_suffix") <- attr(newx, "li_var_suffix")
    x_late <- cbind(
        predict(object = object$model1, newx = x_early, 
            s = object$best_hyperparams$lambda_index),
        newx[, !early_bool]
    )
    y <- predict(object = object$model2, newx = x_late)
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
