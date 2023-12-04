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
  data_spec <- DataSpec()
  model_spec <- ModelSpec(
    fitter = zeroSum::zeroSum,
    response_type = "survival_censored",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    pfs_leq = 1.8
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
  model_spec$pfs_leq <- 2.3
  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
    )
  )
})
