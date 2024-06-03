test_that("nested_pseudo_cv() works", {

    set.seed(123)

    n_samples <- 100
    n_genes <- 10
    n_fold <- 3
    lambda <- c(1, 2)

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y
    hyperparams1 <- list(family = "binomial", nFold = n_fold, lambda = lambda, 
        zeroSum = FALSE)
    hyperparams2 <- list(mtry = c(3, 4, ncol(x)+1), min.node.size = c(4, 5), 
        classification = TRUE, num.trees = 100, skip_on_invalid_input = TRUE)

    fit <- nested_pseudo_cv(
        x = x, y = y, fitter1 = ptk_zerosum, fitter2 = ptk_ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        metric = "accuracy"
    )
    expect_s3_class(fit, "nested_fit")
    expect_s3_class(fit$model1, "ptk_zerosum")
    expect_s3_class(fit$model2, "ranger")
    hyperparams <- c(hyperparams2, list(lambda = lambda))
    expect_true(all(vapply(
        seq_along(hyperparams),
        function(i) fit$best_hyperparams[[i]] %in% hyperparams[[i]],
        logical(1)
    )))
    expect_equal(dim(fit$search_grid), c(2*(3*2), 2+length(hyperparams2)+1))
    expect_equal(length(fit$model1$val_predict_list), 2)
    expect_true(all(fit$model1$val_predict_list[[1]] >= 0 &
        fit$model1$val_predict_list[[1]] <= 1))

    hyperparams2[["skip_on_invalid_input"]] <- FALSE
    expect_error(nested_pseudo_cv(
        x = x, y = y, fitter1 = ptk_zerosum, fitter2 = ptk_ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        metric = "accuracy"
    ), regexp = "mtry must be less than the number of features.")
    # logistic regression as late model
    # with binomial log-likelihood
    hyperparams2 <- list(family = "binomial", lambda = lambda, zeroSum = FALSE)
    fit <- nested_pseudo_cv(x = x, y = y, fitter1 = ptk_zerosum, fitter2 = 
        ptk_zerosum, hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        metric = "binomial_log_likelihood")
    expect_equal(length(fit$best_hyperparams[["overall_val_performance"]]), 1)
    # with AUC
    x_small <- x[, 1:(n_genes+1)]
    attr(x_small, "li_var_suffix") <- attr(x, "li_var_suffix")
    fit <- nested_pseudo_cv(x = x_small, y = y, fitter1 = ptk_zerosum, fitter2 = 
        ptk_zerosum, hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        metric = "roc_auc")
    metric_v <- fit$search_grid[["overall_model_performance"]]
    expect_true(all(metric_v >= 0 & metric_v <= 1))

    # Errors
    x_small <- x[1:(n_samples-1), ]
    attr(x_small, "li_var_suffix") <- attr(x, "li_var_suffix")
    expect_error(nested_pseudo_cv(
        x = x_small, y = y, fitter1 = ptk_zerosum, fitter2 = ptk_ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2
    ), regexp = "is not TRUE")
    attr(x, "li_var_suffix") <- "--"
    expect_error(nested_pseudo_cv(
        x = x, y = y, fitter1 = ptk_zerosum, fitter2 = ptk_ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2 
    ), regexp = "All features are for the early model")
})

test_that("nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 3
    lambda <- 1
    search_grid <- expand.grid(list(lambda = 1:2, mtry = 3:4))

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y

    fit1 <- structure(1, class = c("zeroSum", "list")) 
    fit2 <- structure(2, class = c("ranger", "list")) 
    fit <- nested_fit(model1 = fit1, model2 = fit2, search_grid = search_grid, 
        best_hyperparams = list(lambda = 123))
    expect_s3_class(fit, "nested_fit")

    # Errors
    class(fit1) <- "nonesense"
    expect_error(nested_fit(model1 = fit1, model2 = fit2,
        search_grid = search_grid, best_hyperparams = list(lambda = lambda, 
        mtry = 3, min.node.size = 4, classification = TRUE, num.trees = 100))
    )
})

test_that("predict.nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 1 
    lambda <- 1
    search_grid <- expand.grid(list(num.trees = 100, mtry = 3, lambda = lambda))

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y
    x_early <- x[, seq(n_genes)]
    attr(x_early, "li_var_suffix") <- attr(x, "li_var_suffix")
    x_late <- x[, -seq(n_genes)]

    fit1 <- ptk_zerosum(x = x_early, y = y, nFold = n_fold, lambda = lambda)
    fit2 <- ptk_ranger(x = cbind(fit1$cv.predict[[1]], x_late), y = y, 
        mtry = 3, min.node.size = 4, classification = TRUE, num.trees = 100) 
    n_fit <- nested_fit(fit1, fit2, search_grid = search_grid, 
        list(lambda_index = seq_along(lambda), lambda = lambda, mtry = 3, 
        min.node.size = 4, classification = TRUE, num.trees = 100)) 
    proj <- predict(n_fit, x)

    expect_equal(length(proj), n_samples)
    expect_true(all(proj %in% c(0, 1)))
    expect_false(is.null(names(proj)))
})

test_that("get_<metric> works", {

    set.seed(432)

    n_samples <- 10

    y <- sample(0:1, n_samples, replace = TRUE)
    y_hat <- sample(0:1, n_samples, replace = TRUE)

    acc <- get_accuracy(y, y_hat)
    expect_true(acc >= 0 & acc <= 1)
    y_hat[1] <- -1.5
    expect_error(get_accuracy(y, y_hat))

    y_hat <- runif(n_samples)
    auc <- get_roc_auc(y, y_hat)
    expect_true(auc >= 0 & auc <= 1)

   bll <- get_binomial_log_likelihood(y, y_hat) 
   y_hat[1] <- 1.1
   expect_error(get_binomial_log_likelihood(y, y_hat))
})
