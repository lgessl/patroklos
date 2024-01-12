test_that("prepare_and_fit", {

  set.seed(4352)

  n_samples <- 30
  n_genes <- 5
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  data_spec <- DataSpec(
    name = "mock",
    directory = "mock_dir",
    train_prop = .5,
    pivot_time_cutoff = 2.,
    cohort = "train"
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]
  pheno_tbl[["split_1"]] <- sample(c("train", "test"), n_samples, replace = TRUE)
  pheno_tbl[["split_2"]] <- sample(c("train", "test"), n_samples, replace = TRUE)
  dir <- withr::local_tempdir()

  model_spec <- ModelSpec(
    name = "cox-zerosum",
    directory = dir,
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, lambda = lambda,
      zeroSum = FALSE),
    response_type = "survival_censored",
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL
  )

  prepare_and_fit(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
  )

  data_spec$cohort <- "test"
  expect_silent(res <- prepare_and_predict(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec,
      lambda = "lambda.min",
      benchmark = "ipi"
  ))

  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(sapply(res, is.list)))
  expect_true(all(sapply(res, length) == length(model_spec$split_index)))
  expect_true(all(sapply(res, function(x) all(sapply(x, is.numeric)))))
})
