test_that("nested_cv() works", {

    set.seed(123)

    n_samples <- 100
    n_genes <- 10
    n_fold <- 3
    lambda <- c(1, 2)
    hyperparams1 <- list(family = "binomial", nfolds = n_fold, lambda = lambda, 
        zeroSum = FALSE)
    hyperparams2 <- list(mtry = c(3, 4), min.node.size = c(4, 5), 
        classification = TRUE, num.trees = 100)
    
    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y

    fit <- nested_cv_oob(
        x = x, y = y, fitter1 = zeroSum::zeroSum, fitter2 = ranger::ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        n_folds = 3, append_to_includes = "++" 
    )
    expect_s3_class(fit, "nested_fit")
    expect_s3_class(fit$model1, "zeroSum")
    expect_s3_class(fit$model2, "ranger")
    hyperparams <- c(list(lambda_index = seq_along(lambda), lambda = lambda), 
        hyperparams2)
    expect_true(all(vapply(
        seq_along(hyperparams),
        function(i) fit$best_hyperparams[[i]] %in% hyperparams[[i]],
        logical(1)
    )))

    # Errors
    x <- x[1:(n_samples-1), ]
    expect_error(nested_cv_oob(
        x = x, y = y, fitter1 = zeroSum::zeroSum, fitter2 = ranger::ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        n_folds = 3, append_to_includes = "++" 
    ))
    append_to_includes <- "--"
    expect_error(nested_cv_oob(
        x = x, y = y, fitter1 = zeroSum::zeroSum, fitter2 = ranger::ranger,
        hyperparams1 = hyperparams1, hyperparams2 = hyperparams2,
        n_folds = 3, append_to_includes = append_to_includes 
    ))
})

test_that("nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 3
    lambda <- 1

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y

    fit1 <- structure(1, class = c("zeroSum", "list")) 
    fit2 <- structure(2, class = c("ranger", "list")) 
    fit <- nested_fit(model1 = fit1, model2 = fit2, best_hyperparams = 
        list(lambda = 123), append_to_includes = "++")
    expect_s3_class(fit, "nested_fit")

    # Errors
    class(fit1) <- "nonesense"
    expect_error(nested_fit(model1 = fit1, model2 = fit2, best_hyperparams = list(lambda = lambda, 
        mtry = 3, min.node.size = 4, classification = TRUE, num.trees = 100), 
        append_to_includes = "--"))
})

test_that("predict.nested_fit() works", {

    set.seed(123)

    n_samples <- 50
    n_genes <- 5
    n_fold <- 1 
    lambda <- 1

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y
    x_early <- x[, seq(n_genes)]
    x_late <- x[, -seq(n_genes)]

    fit1 <- zeroSum::zeroSum(x = x_early, y = y, nFold = n_fold, lambda = lambda)
    fit2 <- ranger::ranger(x = cbind(fit1$cv.predict[[1]], x_late), y = y, 
        mtry = 3, min.node.size = 4, classification = TRUE, num.trees = 100) 
    n_fit <- nested_fit(fit1, fit2, list(lambda_index = seq_along(lambda), 
    lambda = lambda, mtry = 3, min.node.size = 4, classification = TRUE, 
    num.trees = 100), append_to_includes = "++")
    proj <- predict(n_fit, x)

    expect_equal(length(proj), n_samples)
    expect_true(all(proj %in% c(0, 1)))
    expect_false(is.null(names(proj)))
})
