test_that("Prepare works",{

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

  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      model = "cox_lasso_zerosum",
      include_from_continuous_pheno = NULL,
      include_from_discrete_pheno = "discrete_var"
    )
  )
  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      model = "lasso_zerosum",
      include_from_continuous_pheno = "continuous_var",
      include_from_discrete_pheno = "discrete_var",
      pfs_leq = 1.8
    )
  )

  colnames(pheno_tbl)[1] <- "patient"
  colnames(pheno_tbl)[3] <- "pfs"
  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      model = "lasso_zerosum",
      include_from_continuous_pheno = NULL,
      include_from_discrete_pheno = "discrete_var",
      pfs_leq = 2.3,
      patient_id_col = "patient",
      pfs_col = "pfs"
    )
  )
})
