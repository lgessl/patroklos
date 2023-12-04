test_that("prepare_and_fit() works", {

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
  data_spec <- DataSpec()
  dir <- withr::local_tempdir()

  # Case 1: Fit all models specified
  model_spec_1 <- ModelSpec(
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold),
    response_type = "survival_censored",
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    save_dir = file.path(dir, "model1") 
  )
  model_spec_2 <- ModelSpec(
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    save_dir = file.path(dir, "model2"),
    pfs_leq = 2.
  )
  model_spec_list <- list(model_spec_1, model_spec_2)

  expect_message(
    prepare_and_fit(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec_list = model_spec_list
    ),
    regexp = "Creating directory"
  )

  # Case 2: Fit all models specified, but one already exists
  model_spec_3 <- ModelSpec(
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "cox", alpha = 0.5, nFold = n_fold),
    response_type = "survival_censored",
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    save_dir = file.path(dir, "model3") 
  )
  model_spec_list <- list(model_spec_2, model_spec_3)
  expect_message(
    prepare_and_fit(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec_list = model_spec_list
    ),
    regexp = "already exists"
  )
})