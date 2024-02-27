test_that("model$predict() works", {

  set.seed(4352)

  n_samples <- 30
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

  model <- Model$new(
    name = "cox-zerosum",
    directory = dir,
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, lambda = lambda,
      zeroSum = FALSE),
    response_type = "survival_censored"
  )
  model$fit(data, quiet = TRUE)

  data$cohort <- "test"
  res <- model$predict(data, pivot_time_cutoff = 2)

  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(sapply(res, is.list)))
  expect_true(all(sapply(res, length) == length(model$split_index)))
  expect_true(all(sapply(res, function(x) all(sapply(x, is.numeric)))))
})
