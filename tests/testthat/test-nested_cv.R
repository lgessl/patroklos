test_that("nested_pseudo_cv() works", {

    set.seed(123)

    n_samples <- 100
    n_genes <- 10
    n_fold <- 3
    lambda <- c(1, 2)

    xyy <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    hyperparams1 <- list(family = "binomial", nFold = n_fold, lambda = lambda, 
        zeroSum = FALSE)
    hyperparams2 <- list(mtry = c(1, 1.5, n_samples^2), min.node.size = c(4, 5), 
        classification = TRUE, num.trees = 100, skip_on_invalid_input = TRUE)

    fit <- nested_pseudo_cv(
        x = xyy[[1]], y_bin = xyy[[2]], y_cox = xyy[[2]], fitter1 = ptk_zerosum, fitter2 = hypertune(ptk_ranger, 
        metric = "accuracy"), hyperparams1 = hyperparams1, hyperparams2 = hyperparams2
    )
    expect_s3_class(fit, "nested_fit")
    expect_s3_class(fit$model1, "ptk_zerosum")
    expect_s3_class(fit$model2, "ptk_ranger")
    hyperparams <- c(hyperparams2, list(lambda = lambda))
    expect_equal(length(fit$model1$val_predict_list), 2)
    expect_true(all(fit$model1$val_predict_list[[1]] >= 0 &
        fit$model1$val_predict_list[[1]] <= 1))

    # Logistic regression as late model
    hyperparams2 <- list(family = "binomial", lambda = lambda, zeroSum = FALSE)
    fit <- nested_pseudo_cv(x = xyy[[1]], y_bin = xyy[[2]], y_cox = xyy[[3]], 
        fitter1 = ptk_zerosum, fitter2 = ptk_zerosum, hyperparams1 = hyperparams1, 
        hyperparams2 = hyperparams2)

    # Cox and cox
    # y <- cbind(runif(n_samples), sample(c(0, 1), n_samples, replace = TRUE))
    # hyperparams1[["family"]] <- "cox"
    # hyperparams2[["family"]] <- "cox"
    # fit <- nested_pseudo_cv(x = x, y = y, fitter1 = ptk_zerosum, fitter2 = 
    #     ptk_zerosum, hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
    #     metric = "auc")

    # Errors
    attr(xyy[[1]], "li_var_suffix") <- "--"
    expect_error(nested_pseudo_cv(
        x = xyy[[1]], y_bin = xyy[[2]], y_cox = xyy[[2]], fitter1 = ptk_zerosum, 
        fitter2 = ptk_ranger, hyperparams1 = hyperparams1, hyperparams2 = hyperparams2 
    ), regexp = "All features are for the early model")
})

test_that("nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 3
    lambda <- 1
    metric_grid <- matrix(1:9, nrow = 3)
    rownames(metric_grid) <- c("a", "b", "c")
    colnames(metric_grid) <- c("d", "e", "f")

    fit1 <- structure(1, class = c("zeroSum", "list")) 
    fit2 <- structure(2, class = c("ranger", "list")) 
    fit <- nested_fit(model1 = fit1, model2 = fit2, metric_grid = metric_grid, 
        best_hyperparams = list(lambda = 123))
    expect_s3_class(fit, "nested_fit")

    # Errors
    class(fit1) <- "nonesense"
    expect_error(nested_fit(model1 = fit1, model2 = fit2,
        metric_grid = metric_grid, best_hyperparams = list(lambda = lambda, 
        mtry = 3, min.node.size = 4, classification = TRUE, num.trees = 100))
    )
})

test_that("predict.nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 1 
    lambda <- 1

    metric_grid <- matrix(1:9, nrow = 3)
    rownames(metric_grid) <- c("a", "b", "c")
    colnames(metric_grid) <- c("d", "e", "f")
    xyy <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x_early <- xyy[[1]][, seq(n_genes)]
    attr(x_early, "li_var_suffix") <- attr(xyy[[1]], "li_var_suffix")
    x_late <- xyy[[1]][, -seq(n_genes)]

    fit1 <- ptk_zerosum(x = x_early, y_bin = xyy[[2]], y_cox = xyy[[3]], 
        nFold = n_fold, lambda = lambda)
    xx <- intersect_by_names(fit1$val_predict_list[[1]], x_late)
    fit2 <- ptk_ranger(x = cbind(xx[[1]], xx[[2]]), 
        y_bin = xyy[[2]], y_cox = xyy[[3]], mtry = 3, min.node.size = 4, 
        classification = TRUE, num.trees = 100) 
    n_fit <- nested_fit(fit1, fit2, metric_grid = metric_grid, 
        list(lambda_index = seq_along(lambda), lambda = lambda, mtry = 3, 
        min.node.size = 4, classification = TRUE, num.trees = 100)) 
    proj <- predict(n_fit, xyy[[1]])

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
