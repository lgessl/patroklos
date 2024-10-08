#' @title Tune multiple hyperparameters with a single call to a fitter
#' @description Decorate a fitter that can only tune one combination of hyperparameters at a time 
#' into a fitter that can tune multiple combinations at once.
#' @param fitter A fitting function that fulfills the [`fitter_prototype()`] interface.
#' @param select logical. If `TRUE`, the returned `multitune_obj` S3 object will only contain the
#' validated predictions and the parameters of the picked model among the fitted models. 
#' @return A fitter following the interface of [`fitter_prototype()`] documented in 
#' [`multitune_output_prototype()`].
#' @export
multitune <- function(
    fitter,
    select = FALSE) {
    stopifnot(is.logical(select))
    force(fitter)

    cv_fitter <- function(x, y, val_error_fun, ...) {
        # Get all possible combinations of hyperparameters
        grid <- expand.grid(list(...), stringsAsFactors = FALSE)
        fit_obj_list <- lapply(seq(nrow(grid)), function(i) {
            do.call(fitter, c(list(x = x, y = y), as.list(grid[i, ])))
        })
        # Get rid of NAs, i.e. invalid hyperparameters
        na_bool <- is.na(fit_obj_list)
        fit_obj_list <- fit_obj_list[!na_bool]
        grid <- grid[!na_bool, ]
        # Build multitune_obj S3 object
        multitune_obj <- list()
        # Promote validated predictions to top level
        multitune_obj$val_predict_list <- lapply(fit_obj_list, function(x) x$val_predict)
        multitune_obj$lambda <- apply(grid, 1, function(r) {
            paste(names(grid), r,
                sep = "=", collapse = ", "
            )
        })
        # Evaluate according to error
        error_v <- vapply(multitune_obj$val_predict_list, function(y_hat) {
            y_yhat <- intersect_by_names(y[["true"]], y_hat, rm_na = c(TRUE, FALSE))
            val_error_fun(y_yhat[[1]], y_yhat[[2]])
        }, numeric(1))
        multitune_obj$lambda_min_index <- which.min(error_v)
        multitune_obj$lambda_min <- multitune_obj$lambda[multitune_obj$lambda_min_index]
        multitune_obj$val_predict <- multitune_obj$val_predict_list[[multitune_obj$lambda_min_index]]
        multitune_obj$val_error <- error_v
        multitune_obj$min_error <- error_v[multitune_obj$lambda_min_index]
        class(multitune_obj) <- "multitune_obj"
        if (select) {
            fit_obj_list[-multitune_obj$lambda_min_index] <- NA
            multitune_obj$val_predict_list <- NA
        }
        multitune_obj$fit_obj_list <- fit_obj_list
        multitune_obj
    }
}

#' @export
predict.multitune_obj <- function(
    object,
    newx,
    lambda_index = object$lambda_min_index,
    ...) {
    obj <- object$fit_obj_list[[lambda_index]]
    if (length(obj) == 1 && is.na(obj)) {
        stop("Desired model is NA. Did accidentally you set `select = TRUE`?")
    }
    predict(object = obj, newx = newx, ...)
}

# Just there for backward compatibility: renamed ptk_hypertune S3 class into multitune_obj in 
# version 1.0.0
#' @export
predict.ptk_hypertune <- function(
    object,
    newx, 
    lambda_index = object$lambda_min_index,
    ...
){
    class(object) <- "multitune_obj"
    predict(object = object, newx = newx, lambda_index = lambda_index, ...)
}

# A decorator to tune the model-agnostic hyperparameters one can tune for any fitter (or
# for a *uni*versal fitter). Right now, this concerns the `Model` attributes 
# `combine_n_max_categorical_features` and `time_cutoffs`.
unitune <- function(fitter) {
    force(fitter)

    unifitter <- function(x, y_cox, time_cutoffs, combine_n_max_categorical_features,
                          pivot_time_cutoff, val_error_fun, ...) {
        # Get all possible combinations of time_cutoffs, combine_n_max_categorical_features
        grid <- expand.grid(list(
            "time cutoff" = time_cutoffs,
            "n combo max" = combine_n_max_categorical_features
        ), stringsAsFactors = FALSE)

        core <- function(i) {
            n_combo_max <- grid[i, "n combo max"]
            time_cutoff <- grid[i, "time cutoff"]
            x <- trim_combos(x, n_combo_max)
            y <- list()
            y[["bin"]] <- binarize_y(
                y_cox = y_cox, time_cutoff = time_cutoff,
                pivot_time_cutoff = pivot_time_cutoff
            )
            y[["cox"]] <- censor_y(y_cox, time_cutoff)
            y[["true"]] <- binarize_y(
                y_cox = y_cox, time_cutoff =
                    pivot_time_cutoff, pivot_time_cutoff = pivot_time_cutoff
            )
            fit_obj <- fitter(x = x, y = y, val_error_fun, ...)
            fit_obj$combine_n_max_categorical_features <- n_combo_max
            fit_obj$time_cutoff <- time_cutoff
            fit_obj$feature_names <- colnames(x)
            fit_obj
        }

        fit_obj_list <- lapply(seq(nrow(grid)), core)
        errors <- sapply(fit_obj_list, function(x) x$min_error)
        grid[["min error"]] <- errors
        fit_obj <- fit_obj_list[[which.min(errors)]]
        fit_obj$unitune_grid <- grid
        fit_obj
    }
}

#' @title Error rate of a binary classifier
#' @param y A numeric vector with binary entries, the true outcomes.
#' @param y_hat A numeric vector with binary entries, the predicted outcomes.
#' @return Numeric scalar, the error rate.
#' @export
error_rate <- function(y, y_hat) {
    stopifnot(all(y_hat %in% c(0, 1)))
    mean(y != y_hat)
}

#' @title Negative ROC AUC
#' @inheritParams error_rate
#' @param y_hat A numeric vector with continuous entries, the predicted outcomes.
#' @return Numeric scalar, the negative ROC AUC.
#' @export
neg_roc_auc <- function(y, y_hat) {
    pred_obj <- ROCR::prediction(predictions = y_hat, labels = y)
    -ROCR::performance(pred_obj, measure = "auc")@y.values[[1]]
}

#' @title Negative binomial log-likelihood
#' @inheritParams neg_roc_auc
#' @return Numeric scalar, the negative binomial log-likelihood.
#' @export
neg_binomial_log_likelihood <- function(y, y_hat) {
    stopifnot(all(y_hat >= 0 & y_hat <= 1))
    -mean(y * log(y_hat) + (1 - y) * log(1 - y_hat))
}

#' @title Minimal negative precision for thresholds with a minimal prevalence
#' @description A function factory.
#' @param min_prev A numeric scalar, the minimal prevalence.
#' @return A function that takes two numeric vectors `y` and `y_hat` and returns
#' the minimal precision over those thresholds yielding a prevalence of at least
#' `min_prev`.
#' @export
neg_prec_with_prev_greater <- function(min_prev) {
    stopifnot(min_prev >= 0 && min_prev <= 1)
    function(y, y_hat) {
        thresholds <- unique(y_hat)
        prevs <- vapply(thresholds, function(t) mean(y_hat >= t), numeric(1))
        thresholds <- thresholds[prevs >= min_prev]
        precs <- -vapply(thresholds, function(t) mean(y[y_hat >= t]), numeric(1))
        min(precs)
    }
}
