test_that("Prepare() works",{

  set.seed(4352)

  n_samples <- 30
  n_genes <- 5
  n_na_in_pheno <- 5

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]
  data_spec <- DataSpec(
    name = "Mock et al. (2023)",
    directory = "some_dir",
    train_prop = 0.8
  )
  model_spec <- ModelSpec(
    name = "zerosum",
    directory = "some_dir",
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2.,
    response_type = "survival_censored",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var"
  )

  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
    )
  )

  model_spec$response_type <- "binary"
  model_spec$include_from_continuous_pheno <- NULL
  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
    )
  )

  colnames(pheno_tbl)[1] <- "patient"
  colnames(pheno_tbl)[3] <- "pfs"
  data_spec$patient_id_col <- "patient"
  data_spec$pfs_col <- "pfs"
  model_spec$pivot_time_cutoff <- 2.3
  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
    )
  )
})
