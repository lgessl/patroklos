test_that("ptk_ranger() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2
  
  x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")
  x <- x_y$x
  y <- x_y$y

  fit <- ptk_ranger(x = x, y = y, mtry = 1, num.trees = 1, min.node.size = 1, 
    classification = TRUE)
  expect_false(is.null(fit$val_predict))
  expect_true(is.null(fit$predictions))

  expect_error(ptk_ranger(x = x, y = y, mtry = ncol(x)+1, 
    skip_on_invalid_input = FALSE), regexp = "mtry must be")
  
  fit <- ptk_ranger(x = x, y = y, mtry = ncol(x)+1, 
    skip_on_invalid_input = TRUE)
  expect_true(is.na(fit))

  fit <- ptk_ranger(x = x, y = y, mtry = 1.5, rel_mtry = TRUE, 
    num.trees = 1, min.node.size = 1, classification = TRUE)
  expect_equal(fit$mtry, round(sqrt(ncol(x)) * 1.5))  
})

test_that("ptk_zerosum() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2
  lambda <- c(0.5, 1)

  x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")
  x <- x_y$x
  y <- x_y$y

  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = TRUE, 
    binarize_predictions = 0.5, lambda = lambda, nFold = 2, family = "binomial")
  expect_equal(fit$zeroSumWeights, c(rep(1, n_genes), rep(0, 2+1+4)))
  expect_equal(fit$penaltyFactor, c(rep(1, n_genes), rep(0, 2+1+4)))
  expect_equal(fit$binarizePredictions, 0.5)
  expect_false(is.null(fit$val_predict))
  expect_false(is.null(fit$val_predict_list))
  expect_true(is.numeric(fit$val_metric))
  expect_equal(length(fit$val_metric), length(lambda))
  expect_true(is.numeric(fit$best_lambda_index) && length(fit$best_lambda_index) == 1)
  expect_true(is.character(fit$lambda) && length(fit$lambda) == length(lambda))
  expect_true(all(vapply(fit$val_predict_list, function(v) v %in% c(0, 1), 
    logical(n_samples))))
  expect_s3_class(fit, "ptk_zerosum")
  expect_equal(fit$type, 2)

  lambda <- lambda[2]

  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = FALSE, 
    binarize_predictions = NULL, lambda = lambda, nFold = 2, zeroSum = FALSE, 
    family = "binomial")
  expect_false(is.logical(fit$useZeroSum))
  expect_true(is.null(fit$zeroSumWeights))
  expect_true(is.null(fit$penaltyFactor))
  expect_true(is.null(fit$binarizePredictions))
  expect_false(is.null(fit$val_predict))
  expect_true(all(vapply(fit$val_predict_list, function(v) 
    all(v >= 0 & v <= 1), logical(1))), logical(n_samples))
  expect_true(all(fit$val_predict <= 1 & fit$val_predict >= 0))
  expect_equal(fit$type, 2)

  fit <- ptk_zerosum(x = x, y = y, exclude_pheno_from_lasso = TRUE, 
    lambda = lambda, nFold = 2, zeroSum = FALSE, binarize_predictions = 0.6, 
    penalty.factor = runif(ncol(x)), family = "binomial")
  expect_true(all(fit$penaltyFactor[!get_early_bool(x)] == 0))
  expect_true(all(vapply(fit$val_predict_list, function(v) v %in% c(0, 1), 
    logical(n_samples))))
  expect_equal(fit$type, 2)

  expect_error(ptk_zerosum(x, y, binarize_predictions = TRUE))
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
  expect_true(all(y_hat >= 0 & y_hat <= 1))
  
  lambda <- lambda[2]
  y <- cbind(runif(n_samples), sample(c(0, 1), n_samples, replace = TRUE))
  fit <- ptk_zerosum(x = x, y = y, nFold = 2, family = "cox", zeroSum = FALSE)
  y_hat <- predict(fit, x)
  expect_true(all(y_hat > 0))
  expect_equal(fit$type, 4)
})