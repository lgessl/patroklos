#' @title Nest an existing, tuned model together with more features into another 
#' model and tune the latter
#' @description Use the validated predictions of an existing model that only 
#' takes the epression part of the data as input features, and feed them 
#' together with the remaining features into another model. Fit and tune the 
#' second, late model.
#' @inheritParams ptk_zerosum
#' @param x named numeric matrix (samples x features). Features only meant for 
#' the late model are exactly those matching `x`'s `li_var_suffix` attribute.
#' @param model1 `Model` R6 object. The early model trained on the expression 
#' data, with the `fit_obj` attribute set at least in its stored version, i.e., the 
#' early model is already there.
#' @param fitter2 A *patroklos-compliant fitter with CV tuning* (see
#' README for more details).
#' @param hyperparams2 A named list with hyperparameters for the late model.
#' @return A `nested_fit` S3 object.
#' @export
greedy_nestor <- function(
    x,
    y,
    val_error_fun,
    model1,
    fitter2,
    hyperparams2
){
    stopifnot(inherits(model1, "Model"))
    stopifnot(inherits(val_error_fun, "function"))
    stopifnot(inherits(hyperparams2, "list"))

    if (is.null(model1$fit_obj)) {
        model1 <- readRDS(file.path(model1$directory, model1$file))
        if (is.null(model1$fit_obj))
            stop("You need to first train the early model")
    }
    # Backward compatibility
    if (length(class(model1$fit_obj)) == 1 && class(model1$fit_obj)[1] == "list") 
        model1$fit_obj <- model1$fit_obj[[1]]
    
    fit_obj1 <- model1$fit_obj
    x_late <- get_late_x(early_predicted = fit_obj1$val_predict, x = x)
    fit_obj2 <- do.call(fitter2, c(list(x = x_late, y = y, val_error_fun = 
        val_error_fun), hyperparams2))
    error_grid <- as.matrix(fit_obj2$val_error)
    rownames(error_grid) <- fit_obj2$lambda

    nested_fit(model1 = fit_obj1, model2 = fit_obj2, error_grid = error_grid,
        best_hyperparams = list(model1 = model1$lambda_min, model2 = 
        fit_obj2$lambda_min))
}

#' @title Nested cross-validation for second-stage validated predictions
#' @description Perform a nested cross-validation for a late-integration scheme,
#' i.e., perform a cross-validation for the early model and then train second-stage
#' models on validated predictions and evaluate the entire models via 
#' validated predictions made by the late model. "Validated prediction" means 
#' predictions made on independent data like out-of-bag (OOB) or cross-validated 
#' (CV) predictions.
#' @inheritParams ptk_zerosum
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
long_nestor <- function(
    x,
    y,
    val_error_fun,
    fitter1,
    fitter2,
    hyperparams1,
    hyperparams2
){
    # Input checks
    stopifnot(is.matrix(x) && is.numeric(x))
    stopifnot(is.numeric(y[["bin"]]) || is.factor(y[["bin"]]))
    stopifnot(is.numeric(y[["cox"]]) || is.factor(y[["cox"]]))
    if (!all(y[["bin"]] %in% c(0, 1, NA))) {
        stop("y[['bin']] must be binary. But your y[['bin']] has the unique elements: ", 
            paste(unique(y[["bin"]]), collapse = ", "), " (length(y[['bin']]) = ", length(y[["bin"]]), 
            ").")
    }
    stopifnot(is.function(fitter1) && is.function(fitter2))
    stopifnot(is.list(hyperparams1) && is.list(hyperparams2))
    
    x_early <- get_early_x(x = x)

    # First stage
    cv1 <- do.call(
        fitter1,
        c(list(x = x_early, y = y, val_error_fun = function(y, y_hat) 0), hyperparams1)
    )
    # Second stage
    second_cv <- function(i) {
        x_late <- get_late_x(early_predicted = cv1$val_predict_list[[i]], x = x)    
        if (ncol(x_late) != ncol(x)-ncol(x_early)+1)
            stop("Something went wrong with adding the early model's predictions.")
        do.call(fitter2, c(list(x = x_late, y = y, val_error_fun = val_error_fun), 
            hyperparams2))
    }
    cv2_list <- lapply(seq_along(cv1$val_predict_list), second_cv)
    # Choose the best combination
    cv1$lambda_min_index <- which.min(vapply(cv2_list, function(x) x$min_error, 
        numeric(1)))
    cv1$lambda_min <- cv1$lambda[cv1$lambda_min_index]
    cv2 <- cv2_list[[cv1$lambda_min_index]]
    model1 <- cv1
    model2 <- cv2
    if (inherits(cv1, "ptk_hypertune"))
        model1 <- cv1$fit_obj_list[[cv1$lambda_min_index]]
    if (inherits(cv2, "ptk_hypertune")) {
        model2 <- cv2$fit_obj_list[[cv2$lambda_min_index]]
        model2$min_error <- cv2$min_error
    }
    # Report on hypertuning
    # lambda sequence in cv2 may depend on lambda sequence in cv1, in particular 
    # lengths may differ (e.g. because of early stopping), hence we need 
    # to pad the error_grid with -Inf
    cv2_n_lambda_max <- max(sapply(cv2_list, function(x) length(x$lambda)))
    error_grid <- sapply(cv2_list, function(x) {
        c(x$val_error, rep(Inf, cv2_n_lambda_max - length(x$val_error)))
    })
    dim(error_grid) <- c(cv2_n_lambda_max, length(cv1$lambda))
    error_grid <- t(error_grid)
    rownames(error_grid) <- cv1$lambda
    colnames(error_grid) <- paste0("lambda2_", seq(cv2_n_lambda_max)) 

    nested_fit(model1 = model1, model2 = model2, 
        error_grid = error_grid, 
        best_hyperparams = list(model1 = cv1$lambda_min, model2 = cv2$lambda_min))
}

#' @title Construct a nested_fit S3 object
#' @description Construct a nested_fit object that holds the early and late model
#' plus their hyperparameters.
#' @param model1 An S3 object with a `predict` method. The early model.
#' @param model2 An S3 object with a `predict` method. The late model.
#' @param error_grid named numeric matrix. Entry (i, j) is the validated metric 
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
    error_grid,
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
    stopifnot(is.matrix(error_grid))
    stopifnot(!is.null(dimnames(error_grid)))
    stopifnot(is.numeric(error_grid))
    stopifnot(is.list(best_hyperparams))
    structure(
        list(
            "model1" = model1,
            "model2" = model2,
            "val_predict" = model2$val_predict,
            "error_grid" = error_grid,
            "best_hyperparams" = best_hyperparams,
            "min_error" = model2$min_error
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
    rownames(early_predicted) <- rownames(x_early)
    x_late <- get_late_x(early_predicted = early_predicted, x = newx)
    y <- predict(object = object$model2, newx = x_late)
    if (!is.numeric(y) || length(y) != nrow(newx))
        stop("The predict method of `model2` must return a numeric vector of same 
            length as the number of rows in `newx` or a list with the first element 
            being the former.")
    y
}

#' @rdname non_zero_coefs
#' @details Refers to the late model.
#' @export
non_zero_coefs.nested_fit <- function(fit_obj, quiet) {
    non_zero_coefs(fit_obj$model2, quiet = quiet)
}