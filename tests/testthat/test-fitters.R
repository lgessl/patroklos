test_that("ptk_ranger() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2
  
  x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")

  fit <- ptk_ranger(x3y[[1]], x3y[[2]], x3y[[3]], rel_mtry = FALSE, mtry = 1, 
    num.trees = 1, min.node.size = 1, classification = TRUE)
  expect_false(is.null(fit$val_predict))
  expect_true(is.null(fit$predictions))

  expect_error(ptk_ranger(x3y[[1]], x3y[[2]], x3y[[3]], rel_mtry = FALSE, 
    mtry = ncol(x3y[[1]])+1, skip_on_invalid_input = FALSE), regexp = "mtry must be")
  
  fit <- ptk_ranger(x3y[[1]], x3y[[2]], x3y[[3]], rel_mtry = FALSE, 
    mtry = ncol(x3y[[1]])+1, skip_on_invalid_input = TRUE)
  expect_true(is.na(fit))

  fit <- ptk_ranger(x3y[[1]], x3y[[2]], x3y[[3]], rel_mtry = TRUE, mtry = 1.5, 
    num.trees = 5, min.node.size = 1, classification = TRUE)
  expect_equal(fit$mtry, round(sqrt(ncol(x3y[[1]])) * 1.5))  
  y_hat <- predict(fit, x3y[[1]])
  expect_false(all(y_hat %in% c(0, 1)))
  expect_true(all(y_hat >= 0 & y_hat <= 1))
  expect_equal(nrow(y_hat), nrow(x3y[[1]]))
})

test_that("ptk_zerosum() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2
  lambda <- c(0.5, 1)

  x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")

  fit <- ptk_zerosum(x3y[[1]], x3y[2:4], neg_binomial_log_likelihood, 
    exclude_pheno_from_lasso = TRUE, binarize_predictions = 0.5, lambdaSteps = 2, 
    nFold = 2, family = "binomial")
  expect_equal(fit$zeroSumWeights, c(rep(1, n_genes), rep(0, 2+1+4)))
  expect_equal(fit$penaltyFactor, c(rep(1, n_genes), rep(0, 2+1+4)))
  expect_equal(fit$binarizePredictions, 0.5)
  expect_false(is.null(fit$val_predict))
  expect_false(is.null(fit$val_predict_list))
  expect_true(is.numeric(fit$val_error))
  expect_equal(length(fit$val_error), length(lambda))
  expect_true(is.numeric(fit$lambda_min_index) && length(fit$lambda_min_index) == 1)
  expect_true(is.character(fit$lambda) && length(fit$lambda) == length(lambda))
  expect_s3_class(fit, "ptk_zerosum")
  expect_equal(fit$type, 2)

  lambda <- lambda[2]

  fit <- ptk_zerosum(x3y[[1]], x3y[2:4], neg_roc_auc, exclude_pheno_from_lasso = FALSE, 
    binarize_predictions = NULL, lambdaSteps = 2, nFold = 2, zeroSum = FALSE, 
    family = "cox")
  expect_true(all(fit$val_error <= 0 & fit$val_error >= -1))
  expect_false(is.logical(fit$useZeroSum))
  expect_true(is.null(fit$zeroSumWeights))
  expect_true(is.null(fit$penaltyFactor))
  expect_true(is.null(fit$binarizePredictions))
  expect_false(is.null(fit$val_predict))
  expect_equal(fit$type, 4)

  fit <- ptk_zerosum(x3y[[1]], x3y[2:4], neg_roc_auc, exclude_pheno_from_lasso = TRUE, 
    lambdaSteps = 2, nFold = 2, zeroSum = FALSE, binarize_predictions = 0.6, 
    penalty.factor = runif(ncol(x3y[[1]])), family = "binomial")
  expect_true(all(fit$penaltyFactor[!get_early_bool(x3y[[1]])] == 0))
  expect_true(all(vapply(fit$val_predict_list, function(v) v %in% c(0, 1), 
    logical(sum(!is.na(x3y[[2]]))))))
  expect_equal(fit$type, 2)

  expect_error(ptk_zerosum(x3y[[1]], x3y[2:4], neg_roc_auc, binarize_predictions = TRUE))
})

test_that("predict.ptk_zerosum() works", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2

  x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")
  fit <- ptk_zerosum(x3y[[1]], x3y[2:4], neg_binomial_log_likelihood, 
    exclude_pheno_from_lasso = TRUE, binarize_predictions = 0.5, 
    lambdaSteps = 10, nFold = 2, zeroSum = FALSE)

  y_hat <- predict(fit, x3y[[1]], s = 1)
  expect_true(all(y_hat %in% c(0, 1)))

  fit$binarizePredictions <- NULL
  y_hat <- predict(fit, x3y[[1]], s = 2)
  expect_true(all(y_hat >= 0 & y_hat <= 1))

  # Check if predict reacts properly to lambda_min_index attribute
  fit$lambda_min_index <- 2
  y_hat11 <- predict(fit, x3y[[1]])
  fit$lamba_min_index <- 4
  y_hat21 <- predict(fit, x3y[[1]])
  class(fit) <- "zeroSum"
  y_hat12 <- predict(fit, x3y[[1]], s = 2, type = "response")
  y_hat22 <- predict(fit, x3y[[1]], s = 4, type = "response")
  fit$lambdaMinIndex <- 2
  y_hat13 <- predict(fit, x3y[[1]], type = "response")
  fit$lambdaMinIndex <- 4
  y_hat23 <- predict(fit, x3y[[1]], type = "response")
  expect_equal(y_hat11, y_hat12)
  expect_equal(y_hat21, y_hat22)
  expect_equal(y_hat11, y_hat13)
  expect_equal(y_hat21, y_hat23)
  
  y <- cbind(runif(n_samples), sample(c(0, 1), n_samples, replace = TRUE))
  fit <- ptk_zerosum(x3y[[1]], x3y[2:4], neg_roc_auc, nFold = 2, family = "cox", 
    zeroSum = FALSE, lambdaSteps = 10)
  y_hat <- predict(fit, x3y[[1]])
  expect_true(all(y_hat > 0))
  expect_equal(fit$type, 4)
})