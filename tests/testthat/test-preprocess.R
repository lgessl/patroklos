test_that("discretize_ipi_features() works", {

  n_samples <- 5
  pheno_tbl <- tibble::tibble(
    "age" = rnorm(n_samples, 60, 10),
    "ann_arbor_stage" = sample(1:4, n_samples, replace = TRUE),
    "ldh_ratio" = runif(n_samples, 0, 2),
    "ecog_performance_status" = sample(0:4, n_samples, replace = TRUE),
    "n_extranodal_sites" = sample(0:6, n_samples, replace = TRUE)
  )
  data_spec <- DataSpec(name = "test")

  pheno_tbl <- discretize_ipi_features(pheno_tbl, data_spec)
  
  expect_equal(ncol(pheno_tbl), 2*n_samples)
  new_col_bool <- stringr::str_detect(colnames(pheno_tbl), ">")
  expect_equal(sum(new_col_bool), n_samples)
  expect_true(all(as.matrix(pheno_tbl[, new_col_bool]) %in% c(0., 1.)))
})
