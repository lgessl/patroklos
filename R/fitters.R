#' @title Wrap [`ranger::ranger()`] into a patroklos-compliant fitter 
#' @description This function is a patroklos-compliant fitter with validated 
#' predictions, it has a return value with a `val_predict` attribute.
#' @inheritParams fitter_prototype
#' @param mtry integer. Number of features to consider at each split when growing a tree of the 
#' random forest.
#' @param rel_mtry logical. If `TRUE`, interpret `mtry` as relative to
#' `sqrt(ncol(x))` (the recommended value), rounded to the next integer. 
#' Otherwise, take `mtry` directly. 
#' @param skip_on_invalid_input Logical. If `TRUE` and invalid input is detected,
#' return `NA` instead of an error. This is useful when calling this function 
#' from inside [`long_nestor()`].
#' @param ... Further arguments passed to the wrapped function.
#' @return A `ptk_ranger` S3 object, a `ranger` S3 object with the (OOB) 
#' `predictions` attribute renamed to `val_predict`.
#' @export
ptk_ranger <- function(x, y, rel_mtry, mtry = NULL, 
    skip_on_invalid_input = FALSE, ...){
    stopifnot(is.logical(rel_mtry))
    y_bin <- y[["bin"]]

    if (!is.null(mtry)) {
        if (rel_mtry) mtry <- round(sqrt(ncol(x)) * mtry)
        if (ncol(x) < mtry) {
            if (skip_on_invalid_input) return(NA)
            stop("mtry must be less than the number of features.")
        }
    }
    x_y <- intersect_by_names(x, y_bin, rm_na = c(TRUE, TRUE))
    ptk_ranger_obj <- ranger::ranger(x = x_y[[1]], y = x_y[[2]], mtry = mtry, ...)
    # Rename OOB predictions
    val_predict <- as.matrix(ptk_ranger_obj$predictions)
    rownames(val_predict) <- rownames(x_y[[1]])
    ptk_ranger_obj$val_predict <- val_predict
    ptk_ranger_obj$min_error <- ptk_ranger_obj$prediction.error
    ptk_ranger_obj$predictions <- NULL
    ptk_ranger_obj$prediction.error <- NULL
    class(ptk_ranger_obj) <- c("ptk_ranger", class(ptk_ranger_obj))
    ptk_ranger_obj
}

#' @title Wrap [`ranger::predict.ranger()`] into a patroklos-compliant predict
#' function
#' @param object A ptk_ranger S3 object.
#' @param newx Predictor data as a named numeric matrix. Samples correspond to 
#' rows, an attribute `li_var_suffix` is expected to be present.
#' @param ... Further arguments passed to the wrapped function.
#' @return A named numeric vector of ensemble vote share for every sample (
#' i.e. a continuous value between 0 and 1).
#' @export
predict.ptk_ranger <- function(object, newx, ...){
    class(object) <- "ranger"
    y <- predict(object, data = newx, predict.all = TRUE, ...)[[1]]
    as.matrix(rowMeans(y))
}

#' @title Wrap [`zeroSum::zeroSum()`] into a patroklos-compliant fit function
#' @description This function is a patroklos-compliant fitter with integrated 
#' CV and, if `length(lambda) == 1`, also a patroklos-compliant fitter with 
#' validated predictions.
#' @inheritParams zeroSum::zeroSum
#' @param y Named list with entries
#' * `"bin"`, a named numeric one-column matrix, binary response to be used 
#' for training,
#' * `"cox"`, a named numeric two-column matrix to be used for training, time to 
#' event and event (0 = censoring, 1 = event) in first and second column, respectively.
#' * `"true"`, a named numeric one-column matrix, binary response to be used for
#' calculating the CV error.
#' @param val_error_fun Function used to calculate the error of independently 
#' validated predictions. Must take two numeric vector of equal length:
#' `y` and `y_hat`, the true and predicted outcomes, respectively, and 
#' return a numeric scalar; the lower, the better the model. See
#' [`error_rate()`] or [`neg_roc_auc()`] for examples. 
#' @param exclude_pheno_from_lasso Logical. If `TRUE`, set LASSO penalty weights
#' corresponding to features from the pheno data to zero.
#' @param binarize_predictions numeric or NULL. If not NULL, the predict method 
#' for the returned `ptk_zerosum` object will binarize the predictions using the
#' `binarize_predictions` as a threshold.
#' @return A `ptk_zerosum` S3 object with the `cv_predict` attribute renamed to 
#' `val_predict_list`.
#' @export
ptk_zerosum <- function(
    x,
    y,
    val_error_fun,
    exclude_pheno_from_lasso = TRUE,
    binarize_predictions = NULL,    
    ...,
    nFold = 10,
    zeroSum.weights = NULL, # Do not partially match these
    penalty.factor = NULL,
    family = "binomial"
){
    stopifnot(is.logical(exclude_pheno_from_lasso))
    stopifnot(is.null(binarize_predictions) || 
        (is.numeric(binarize_predictions) && binarize_predictions > 0 && 
        binarize_predictions < 1))

    early_bool <- get_early_bool(x, for_li = FALSE)
    # zeroSum weights
    if (is.null(zeroSum.weights))
        zeroSum.weights <- as.numeric(early_bool) 
    # LASSO penalty factors
    if (exclude_pheno_from_lasso) {
        if(is.null(penalty.factor)) {
            penalty.factor <- as.numeric(early_bool)
        } else {
            penalty.factor[!early_bool] <- 0
        }
    }

    y4f <- y[["bin"]] # y for fitter
    # Finally prepare y
    if (family == "cox") y4f <- y[["cox"]]
    x_y <- intersect_by_names(x, y4f, rm_na = c(TRUE, TRUE))
    x <- x_y[[1]]
    y4f <- x_y[[2]]
    if (family == "cox") {
        # zeroSum removes censored samples with time lower than the first event
        ord <- order(y4f[, 1], y4f[, 2])
        i <- 1
        while (y4f[ord[i], 2] == 0) i <- i + 1
        ord <- ord[i:length(ord)]
        x <- x[ord, , drop = FALSE]
        y4f <- y4f[ord, , drop = FALSE]
    } 

    fit_obj <- zeroSum::zeroSum(x = x, y = y4f, nFold = nFold, 
        zeroSum.weights = zeroSum.weights, penalty.factor = penalty.factor, 
        family = family, ...)

    fit_obj$val_predict_list <- lapply(fit_obj$cv_predict, function(y_hat) {
        rownames(y_hat) <- rownames(y4f)
        if (family == "binomial") y_hat <- 1/(1+exp(-y_hat))
        y_hat
    })
    fit_obj$cv_predict <- NULL
    fit_obj$val_error <- vapply(fit_obj$val_predict_list, function(y_hat) {
        y_yhat <- intersect_by_names(y[["true"]], y_hat, rm_na = c(TRUE, FALSE))
        val_error_fun(y_yhat[[1]], y_yhat[[2]])
    }, numeric(1))
    fit_obj$lambda_min_index <- which.min(fit_obj$val_error)
    fit_obj$lambdaMinIndex <- fit_obj$lambda_min_index
    # Lambda seq is often too long because of early stopping
    fit_obj$lambda <- fit_obj$lambda[seq_along(fit_obj$val_predict_list)]
    fit_obj$lambda <- paste("lambda", fit_obj$lambda, sep = "=")

    if (is.numeric(binarize_predictions))
        fit_obj$val_predict_list <- lapply(fit_obj$val_predict_list,
            function(v) 1 * (v > binarize_predictions))
    fit_obj$val_predict <- fit_obj$val_predict_list[[fit_obj$lambda_min_index]]
    if (fit_obj$useZeroSum)
        fit_obj$zeroSumWeights <- zeroSum.weights
    if (fit_obj$standardize == 0)
        fit_obj$penaltyFactor <- penalty.factor
    fit_obj$binarizePredictions <- binarize_predictions
    fit_obj$min_error <- fit_obj$val_error[fit_obj$lambda_min_index]
    class(fit_obj) <- "ptk_zerosum"
    return(fit_obj)
}

#' @title Wrap [`zeroSum::predict.zeroSum()`] into a patroklos-compliant predict function
#' @description `Model$predict()` uses this function as a predict method for
#' the `predict()` generic.
#' @inheritParams predict.ptk_ranger
#' @param object A ptk_patroklos S3 object.
#' @param newx Predictor data as a *named* numeric matrix. Samples correspond to 
#' rows.
#' @return A numeric vector with the prediction for every sample. Unlike 
#' `zeroSum::predict.zeroSum()`, 
#' * if `object$type == 2` (i.e., `object` was fitted with `family = "binomial"
#' by `ptk_zerosum`), the predictions are binomial probabilities,
#' * if `object$type == 4` (i.e., `object` was fitted with `family = "cox" by
#' `ptk_zerosum`), the predictions are the exponential of the linear predictor, 
#' i.e., `exp(beta_0 + beta^T x)`.
#' * If `object$binarizePredictions` is not `NULL`, the predictions are binarized
#' via `as.numeric(y > object$binarizePredictions)`.
#' @export
predict.ptk_zerosum <- function(
    object,
    newx,
    ...
){
    class(object) <- "zeroSum"
    args <- list(...)
    if(!is.null(args[["type"]]))
        stop("`type` argument is not allowed.") 
    stopifnot(all(object$variables.names[-1] == colnames(newx)))
    type <- NULL
    if(object$type %in% c(2, 4)) # 2 \mapsto binomial, 4 \mapsto cox
        # get probabilities (binomial), exp(beta_0 + beta^T x) (cox)
        type <- "response"     
    object$lambdaMinIndex <- object$lambda_min_index
    y <- predict(object = object, newx = newx, type = type, ...)
    if (!is.null(object$binarizePredictions))
        y <- 1 * (y > object$binarizePredictions)
    return(y)
}

#' @title Project on a feature
#' @description Trivially wraps `ptk_zerosum` to project on a feature, usually the output of a 
#' benchmark model.
#' @inheritParams ptk_zerosum
#' @param feature character. Project on this feature.
#' @return A `ptk_zerosum` S3 object with the `val_predict` attribute set to the predictions on 
#' the "train" cohort. The benchmark usually has not been trained on this train cohort such that 
#' `val_predict` are predictions on independent data.
#' @export
projection_on_feature <- function(
    x, 
    y, 
    val_error_fun,
    feature,
    ...
){
    feature <- paste0(feature, attr(x, "li_var_suffix"))
    if (!feature %in% colnames(x)) stop("Feature not found in predictor matrix")
    y[["bin"]] <- x[, feature, drop = FALSE]
    fit_obj <- ptk_zerosum(x = x, y = y, val_error_fun = val_error_fun, lambda = 0, nFold = 1, 
        family = "gaussian")
    fit_obj$coef[[1]][, 1] <- 1 * c(0, colnames(x) == feature)
    fit_obj$val_predict <- y[["bin"]]
    return(fit_obj)
}

#' @title Get features with non-zero coefficients
#' @param fit_obj S3 fit object.
#' @return One-column matrix with non-zero coefficients and featue names as row names.
#' @details Throws an error if method for a class is not implemented (e.g. for `ptk_ranger` 
#' because "non-zero coefficients" do not make sense with random forests).
#' @export
non_zero_coefs <- function(fit_obj){
    UseMethod("non_zero_coefs")
}

#' @rdname non_zero_coefs
#' @export
non_zero_coefs.ptk_zerosum <- function(fit_obj){
    coefs <- fit_obj$coef[[fit_obj$lambda_min_index]]
    rownames(coefs) <- fit_obj$variables.names
    m <- coefs[coefs[, 1] != 0, , drop = FALSE]
}