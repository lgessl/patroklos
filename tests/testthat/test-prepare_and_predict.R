test_that("prepare_and_fit", {

  set.seed(4352)

  n_samples <- 30
  n_genes <- 5
  n_na_in_pheno <- 5
  n_fold <- 2

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]
  data_spec <- DataSpec(name = "Mock et al. (2023)")
  dir <- withr::local_tempdir()

  model_spec_1 <- ModelSpec(
    name = "cox-zerosum",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold),
    response_type = "survival_censored",
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    save_dir = file.path(dir, "cox"),
    pfs_leq = 2.
  )
  model_spec_2 <- ModelSpec(
    name = "binomial-zerosum",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    save_dir = file.path(dir, "binomial"),
    pfs_leq = 2.
  )
  model_spec_list <- list(model_spec_1, model_spec_2)

  prepare_and_fit(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec_list = model_spec_list
  )

  expect_silent(pred_obs <- prepare_and_predict(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec_1
  ))
  expect_true(is.list(pred_obs))
  expect_true(is.vector(pred_obs[["predicted"]]))
  expect_true(is.vector(pred_obs[["actual"]]))
  expect_true(is.numeric(pred_obs[["predicted"]]))
  expect_true(is.numeric(pred_obs[["actual"]]))

  expect_silent(pred_obs <- prepare_and_predict(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec_2
  ))
})
