test_that("model$predict() works", {

  set.seed(4352)

  n_samples <- 50
  n_genes <- 5
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )

  model1 <- Model$new(
    name = "cox-zerosum",
    directory = dir,
    fitter = ptk_zerosum,
    split_index = 1:2,
    time_cutoffs = Inf,
    val_error_fun = neg_roc_auc,
    hyperparams = list(family = "cox", alpha = 1, nFold = n_fold, 
    lambda = lambda, zeroSum = FALSE)
  )
  model1$fit(data, quiet = TRUE)
  model2 <- Model$new(
    name = "logistic",
    directory = dir,
    fitter = ptk_zerosum,
    split_index = 1:2,
    time_cutoffs = 2,
    val_error_fun = neg_binomial_log_likelihood,
    hyperparams = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE)
  )
  model2$fit(data, quiet = TRUE)

  data$cohort <- "test"
  res <- model1$predict(data, quiet = TRUE)
  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(sapply(res, is.list)))
  expect_true(all(sapply(res, length) == length(model1$split_index)))
  expect_true(all(sapply(res, function(x) all(sapply(x, is.numeric)))))
  expect_equal(length(table(res[["actual"]][[2]])), 2)

  data$cohort <- "val_predict"
  res <- model2$predict(data, quiet = TRUE)
  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(sapply(res[1:2], is.list))) # third list element is NULL
  expect_true(all(sapply(res[1:2], length) == length(model2$split_index)))
  expect_true(all(sapply(res, function(x) all(sapply(x, is.numeric)))))
  expect_equal(length(table(res[["actual"]][[2]])), 2)
})
