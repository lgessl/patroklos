test_that("greedy_nestor() works", {

    set.seed(123)

    n_samples <- 100
    n_genes <- 5
    n_fold <- 2
    n_lambda <- 2
    dir <- withr::local_tempdir()

    x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    data <- generate_mock_data(n_samples = n_samples, n_genes = n_genes)
    hyperparams1 <- list(family = "binomial", nFold = n_fold, lambdaStep = n_lambda,
        zeroSum = FALSE)
    hyperparams2 <- list(rel_mtry = TRUE, mtry = c(1, 1.5, n_samples^2), 
        min.node.size = c(4, 5), classification = TRUE, num.trees = 100, 
        skip_on_invalid_input = TRUE)
    
    model1 <- Model$new(name = "log", fitter = ptk_zerosum, hyperparams = hyperparams1,
        time_cutoffs = 2, val_error_fun = neg_prec_with_prev_greater(0.15), 
        directory = dir)
    model1$fit(data, quiet = TRUE)
    model1$fit_obj <- NULL
    nested_fit <- greedy_nestor(x = x3y[[1]], y = x3y[2:4], val_error_fun = 
        neg_prec_with_prev_greater(0.15), model1 = model1, fitter2 = hypertune(ptk_ranger),
        hyperparams2 = hyperparams2)
    yhat <- predict(nested_fit, x3y[[1]])
    expect_s3_class(nested_fit, "nested_fit")
    expect_s3_class(nested_fit$model1, "ptk_zerosum")
    expect_s3_class(nested_fit$model2, "ptk_hypertune")
    expect_true(all(nested_fit$val_predict %in% c(0, 1)))
    expect_true(all(yhat <= 1 & yhat >= 0))

    hyperparams2 <- list(family = "binomial", nFold = n_fold, lambdaStep = n_lambda,
        zeroSum = FALSE)
    nested_fit <- greedy_nestor(x = x3y[[1]], y = x3y[2:4], val_error_fun = 
        neg_prec_with_prev_greater(0.15), model1 = model1, fitter2 = ptk_zerosum,
        hyperparams2 = hyperparams2)
    yhat <- predict(nested_fit, x3y[[1]])
    expect_s3_class(nested_fit, "nested_fit")
    expect_s3_class(nested_fit$model1, "ptk_zerosum")
    expect_s3_class(nested_fit$model2, "ptk_zerosum")
    expect_true(all(nested_fit$val_predict <= 1 & nested_fit$val_predict >= 0))
})

test_that("long_nestor() works", {

    set.seed(123)

    n_samples <- 100
    n_genes <- 10
    n_fold <- 2
    lambda <- c(1, 2)

    x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    hyperparams1 <- list(family = "binomial", nFold = n_fold, lambda = lambda, 
        zeroSum = FALSE)
    hyperparams2 <- list(rel_mtry = TRUE, mtry = c(1, 1.5, n_samples^2), 
        min.node.size = c(4, 5), classification = TRUE, num.trees = 100, 
        skip_on_invalid_input = TRUE)

    fit <- long_nestor(
        x = x3y[[1]], x3y[2:4], val_error_fun = error_rate, fitter1 = ptk_zerosum, 
        fitter2 = hypertune(ptk_ranger), hyperparams1 = 
        hyperparams1, hyperparams2 = hyperparams2
    )
    expect_s3_class(fit, "nested_fit")
    expect_s3_class(fit$model1, "ptk_zerosum")
    expect_s3_class(fit$model2, "ptk_ranger")
    hyperparams <- c(hyperparams2, list(lambda = lambda))
    expect_equal(length(fit$model1$val_predict_list), 2)
    expect_true(all(fit$model1$val_predict_list[[1]] >= 0 &
        fit$model1$val_predict_list[[1]] <= 1))

    # Logistic regression as late model
    hyperparams1 <- list(family = "binomial", lambdaSteps = 3, zeroSum = FALSE, 
        nFold = n_fold)
    hyperparams2 <- list(family = "binomial", lambdaSteps = 3, zeroSum = FALSE, 
        nFold = n_fold)
    fit <- long_nestor(x = x3y[[1]], x3y[2:4], val_error_fun = neg_binomial_log_likelihood,
        fitter1 = ptk_zerosum, fitter2 = ptk_zerosum, hyperparams1 = hyperparams1, 
        hyperparams2 = hyperparams2)

    # Cox and cox
    hyperparams1 <- list(family = "cox", lambdaSteps = 2, zeroSum = FALSE,
        nFold = n_fold)
    hyperparams2 <- list(family = "cox", lambdaSteps = 2, zeroSum = FALSE, 
        nFold = n_fold)
    fit <- long_nestor(x = x3y[[1]], y = x3y[2:4], val_error_fun = neg_roc_auc,
        fitter1 = ptk_zerosum, fitter2 = ptk_zerosum, hyperparams1 = hyperparams1, 
        hyperparams2 = hyperparams2)

    # Binarized logistic and rf
    hyperparams1 <- list(family = "binomial", lambdaSteps = 2, zeroSum = FALSE, 
        nFold = n_fold, binarize_predictions = 0.5)
    hyperparams2 <- list(rel_mtry = FALSE, mtry = 1, min.node.size = c(4, 5), 
        classification = TRUE, num.trees = 100, skip_on_invalid_input = TRUE)
    fit <- long_nestor(x = x3y[[1]], y = x3y[2:4], val_error_fun = error_rate,
        fitter1 = ptk_zerosum, fitter2 = hypertune(ptk_ranger), 
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2)

    # Errors
    attr(x3y[[1]], "li_var_suffix") <- "--"
    expect_error(long_nestor(
        x = x3y[[1]], y = x3y[2:4], val_error_fun = neg_roc_auc, fitter1 = ptk_zerosum, 
        fitter2 = ptk_ranger, hyperparams1 = hyperparams1, hyperparams2 = hyperparams2 
    ), regexp = "All features are for the early model")
})

test_that("nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 3
    lambda <- 1
    error_grid <- matrix(1:9, nrow = 3)
    rownames(error_grid) <- c("a", "b", "c")
    colnames(error_grid) <- c("d", "e", "f")

    fit1 <- structure(1, class = c("zeroSum", "list")) 
    fit2 <- structure(list(min_error = 1.1), class = c("ranger", "list")) 
    fit <- nested_fit(model1 = fit1, model2 = fit2, error_grid = error_grid, 
        best_hyperparams = list(lambda = 123))
    expect_s3_class(fit, "nested_fit")

    # Errors
    class(fit1) <- "nonesense"
    expect_error(nested_fit(model1 = fit1, model2 = fit2,
        error_grid = error_grid, best_hyperparams = list(lambda = lambda, 
        mtry = 3, min.node.size = 4, classification = TRUE, num.trees = 100))
    )
})

test_that("predict.nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 1 

    error_grid <- matrix(1:9, nrow = 3)
    rownames(error_grid) <- c("a", "b", "c")
    colnames(error_grid) <- c("d", "e", "f")
    x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x_early <- x3y[[1]][, seq(n_genes)]
    attr(x_early, "li_var_suffix") <- attr(x3y[[1]], "li_var_suffix")
    x_late <- x3y[[1]][, -seq(n_genes)]

    fit1 <- ptk_zerosum(x = x_early, y = x3y[2:4], val_error_fun = neg_roc_auc,
        nFold = n_fold, lambdaSteps = 5)
    xx <- intersect_by_names(fit1$val_predict_list[[1]], x_late)
    fit2 <- ptk_ranger(x = cbind(xx[[1]], xx[[2]]), 
        rel_mtry = FALSE, y_bin = x3y[[2]], y_cox = x3y[[3]], mtry = 3, 
        min.node.size = 4, classification = TRUE, num.trees = 100) 
    n_fit <- nested_fit(fit1, fit2, error_grid = error_grid, best_hyperparams = 
        list(model1 = "lambda=1.5", model2 = "mtry = 3"))
    proj <- predict(n_fit, x3y[[1]])

    expect_equal(length(proj), n_samples)
    expect_true(all(proj <= 1 & proj >= 0))
    expect_false(all(proj %in% c(0, 1)))

    # Check if predictions are sensitive to lambda_min_index
    x_late <- cbind(xx[[1]], xx[[2]])
    attr(x_late, "li_var_suffix") <- attr(x3y[[1]], "li_var_suffix")
    fit2 <- ptk_zerosum(x = x_late, y = x3y[2:4], val_error_fun = neg_roc_auc,
        nFold = n_fold, lambdaSteps = 2)
    # Make sure second model actually uses prediction of model 1, and that it 
    # does so differently
    fit2$coef[[1]][2, 1] <- 1
    fit2$coef[[2]][2, 1] <- -1
    n_fit <- nested_fit(fit1, fit2, error_grid = error_grid, best_hyperparams = 
        list(model1 = "lambda=1.5", model2 = "lambda=1.5"))
    n_fit$model1$lambda_min_index <- 1
    proj1 <- predict(n_fit, x3y[[1]])
    n_fit$model1$lambda_min_index <- 2
    proj2 <- predict(n_fit, x3y[[1]])
    expect_false(all(proj1 == proj2))
})

test_that("error functions work", {

    set.seed(432)

    n_samples <- 10

    y <- sample(0:1, n_samples, replace = TRUE)
    y_hat <- sample(0:1, n_samples, replace = TRUE)

    acc <- error_rate(y, y_hat)
    expect_true(acc >= 0 & acc <= 1)
    y_hat[1] <- -1.5
    expect_error(error_rate(y, y_hat))

    y_hat <- runif(n_samples)
    auc <- neg_roc_auc(y, y_hat)
    expect_true(auc >= -1 & auc <= 0)

   bll <- neg_binomial_log_likelihood(y, y_hat) 
   y_hat[1] <- 1.1
   expect_error(neg_binomial_log_likelihood(y, y_hat))
})
