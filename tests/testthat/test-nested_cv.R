test_that("nested_cv() works", {

    set.seed(123)

    n_samples <- 100
    n_genes <- 10
    n_fold <- 3
    lambda <- c(1, 2)

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
        return_type = "fitter")
    x <- x_y$x
    y <- x_y$y

    fit <- nested_cv_oob(
        x = x, y = y, fitter1 = zeroSum::zeroSum, fitter2 = ranger::ranger,
        hyperparams1 = list(family = "binomial", nfolds = n_fold, lambda = lambda, 
            zeroSum = FALSE),
        hyperparams2 = list(mtry = c(3, 4), min.node.size = c(4, 5), classification = TRUE),
        n_folds = 3, append_to_includes = "++" 
    )
})