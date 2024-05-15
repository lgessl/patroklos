test_that("ptk_zerosum() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2
  lambda <- 1

  x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")
  x <- x_y$x
  y <- x_y$y

  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = TRUE, 
    binarize_predictions = 0.5, lambda = lambda, nFold = 2)
  expect_equal(fit$zeroSumWeights, c(rep(1, n_genes), rep(0, 2+1+4)))
  expect_equal(fit$penaltyFactor, c(rep(1, n_genes), rep(0, 2+1+4)))
  expect_equal(fit$binarizePredictions, 0.5)
  expect_false(is.null(fit$cv_predict))
  expect_false(is.null(fit$cv_predict_list))
  expect_s3_class(fit, "ptk_zerosum")

  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = FALSE, 
    binarize_predictions = NULL, lambda = lambda, nFold = 2, zeroSum = FALSE)
  expect_false(as.logical(fit$useZeroSum))
  expect_true(is.null(fit$zeroSumWeights))
  expect_true(is.null(fit$penaltyFactor))
  expect_true(is.null(fit$binarizePredictions))

  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = TRUE, 
    lambda = lambda, nFold = 2, zeroSum = FALSE, binarize_predictions = 0.6, 
    penalty.factor = runif(ncol(x)))
  expect_true(all(fit$penaltyFactor[!get_early_bool(x)] == 0))
})

test_that("predict.ptk_zerosum() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2
  lambda <- c(0.5, 1)

  x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")
  x <- x_y$x
  y <- x_y$y
  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = TRUE, 
    binarize_predictions = 0.5, lambda = lambda, nFold = 2, zeroSum = FALSE)

  y_hat <- predict(fit, x, s = 1)
  expect_true(all(y_hat %in% c(0, 1)))

  fit$binarizePredictions <- NULL
  y_hat <- predict(fit, x, s = 2)
  expect_false(all(y_hat %in% c(0, 1)))
})