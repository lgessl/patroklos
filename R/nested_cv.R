#' @title Nested cross-validation for second-stage OOB predictions
#' @description Perform a nested cross-validation for a late-integration scheme,
#' i.e., perform a cross-validation for the early model and then train second-stage
#' models on validated predictions and evaluate the entire models via 
#' validated predictions made by the late model. "Validated prediction" means 
#' predictions made on independent data like out-of-bag (OOB) or cross-validated 
#' (CV) predictions.
#' @param x A numeric matrix holding the predictor features: rows are samples and
#' columns are features.
#' @param y A numeric vector holding the response variable. 
#' @param fitter1 A *patroklos-compliant fitter with CV tuning* (see README for 
#' more details).
#' @param fitter2 A *patroklos-compliant fitter with validated predictions* (see 
#' README for more details). If it returns `"next"`, we skip the current 
#' combination of hyperparameters and set its metric to `-Inf`.
#' @param hyperparams1 A named list with hyperparaters we will pass to `fitter1`.
#' @param hyperparams2 A named list with hyperparameters for the late model. 
#' Unlike `hyperparams1`, we call `fitter2` for every combination of values in 
#' `hyperparams2` and lambda value from `fitter1`.
#' @return An S3 object with class `nested_fit`, the model with the best 
#' performance according to validated predictions assessed with `metric`. 
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
    y_bin,
    y_cox,
    fitter1,
    fitter2,
    hyperparams1,
    hyperparams2
){
    # Input checks
    stopifnot(is.matrix(x) && is.numeric(x))
    stopifnot(is.numeric(y_bin) || is.factor(y_bin))
    stopifnot(is.numeric(y_cox) || is.factor(y_cox))
    if (!all(y_bin %in% c(0, 1, NA))) {
        stop("y_bin must be binary. But your y_bin has the unique elements: ", 
            paste(unique(y_bin), collapse = ", "), " (length(y_bin) = ", length(y_bin), 
            ").")
    }
    stopifnot(is.function(fitter1) && is.function(fitter2))
    stopifnot(is.list(hyperparams1) && is.list(hyperparams2))
    
    x_early <- get_early_x(x = x)

    # First stage
    cv1 <- do.call(
        fitter1,
        c(list(x = x_early, y_bin = y_bin, y_cox = y_cox), hyperparams1)
    )
    # Second stage
    second_cv <- function(i) {
        x_late <- get_late_x(early_predicted = cv1$val_predict_list[[i]], x = x)    
        if (ncol(x_late) != ncol(x)-ncol(x_early)+1)
            stop("Something went wrong with adding the early model's predictions.")
        do.call(fitter2, c(list(x = x_late, y_bin = y_bin, y_cox = y_cox), 
            hyperparams2))
    }
    cv2_list <- lapply(seq_along(cv1$val_predict_list), second_cv)
    # Choose the best combination
    cv1$best_lambda_index <- which.max(vapply(cv2_list, function(x) x$best_metric, 
        numeric(1)))
    cv1$best_lambda <- cv1$lambda[cv1$best_lambda_index]
    cv2 <- cv2_list[[cv1$best_lambda_index]]
    model1 <- cv1
    model2 <- cv2
    if (inherits(cv1, "ptk_hypertune"))
        model1 <- cv1$fit_obj_list[[cv1$best_lambda_index]]
    if (inherits(cv2, "ptk_hypertune"))
        model2 <- cv2$fit_obj_list[[cv2$best_lambda_index]]
    # Report on hypertuning
    metric_grid <- sapply(cv2_list, function(x) x$val_metric)
    dim(metric_grid) <- c(length(cv2$lambda), length(cv1$lambda))
    metric_grid <- t(metric_grid)
    rownames(metric_grid) <- cv1$lambda
    colnames(metric_grid) <- cv2$lambda

    nested_fit(model1 = model1, model2 = model2, 
        metric_grid = metric_grid, 
        best_hyperparams = list(model1 = cv1$best_lambda, model2 = cv2$best_lambda))
}

#' @title Construct a nested_fit S3 object
#' @description Construct a nested_fit object that holds the early and late model
#' plus their hyperparameters.
#' @param model1 An S3 object with a `predict` method. The early model.
#' @param model2 An S3 object with a `predict` method. The late model.
#' @param metric_grid named numeric matrix. Entry (i, j) is the validated metric 
#' of the nested model with the early model having the i-th and the late model 
#' having the j-th of their respective hyperparameters (which are given as 
#' dimnames).
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
    metric_grid,
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
    stopifnot(is.matrix(metric_grid))
    stopifnot(!is.null(dimnames(metric_grid)))
    stopifnot(is.numeric(metric_grid))
    stopifnot(is.list(best_hyperparams))
    structure(
        list(
            "model1" = model1,
            "model2" = model2,
            "metric_grid" = metric_grid,
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
    x_early <- get_early_x(x = newx)
    early_predicted <- predict(object = object$model1, newx = x_early)
    x_late <- get_late_x(early_predicted = early_predicted, x = newx)
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

get_accuracy <- function(y, y_hat){
    stopifnot(all(y_hat %in% c(0, 1)))
    mean(y == y_hat)
}

get_roc_auc <- function(y, y_hat){
    pred_obj <- ROCR::prediction(predictions = y_hat, labels = y)
    ROCR::performance(pred_obj, measure = "auc")@y.values[[1]]
}

get_binomial_log_likelihood <- function(y, y_hat){
    stopifnot(all(y_hat > 0) & all(y_hat < 1))
    sum(y * log(y_hat) + (1-y) * log(1-y_hat))
}