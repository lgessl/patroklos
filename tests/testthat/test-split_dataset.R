test_that("split_dataset works", {

  set.seed(234)

  n_samples <- 100
  n_genes <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = 0,
    to_csv = dir
  )
  expr_tbl <- data[["expr_tbl"]]
  pheno_tbl <- data[["pheno_tbl"]]
  data_spec <- DataSpec(
    name = "mock",
    directory = dir
  )

  for(based_on_pfs_cut in c(TRUE, FALSE)){

    split <- split_dataset(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      train_prop = 0.5,
      pfs_cut = 2,
      based_on_pfs_cut = based_on_pfs_cut
    )
    for(part in c("train", "test")){
      expect_no_error(
        qc_preprocess(
          expr_tbl = split[[part]][["expr"]],
          pheno_tbl = split[[part]][["pheno"]],
          data_spec = split[[part]][["data_spec"]]
        )
      )
      expect_equal(
        nrow(split[[part]][["expr"]]),
        nrow(expr_tbl)
      )
      expect_equal(
        colnames(split[[part]][["pheno"]]),
        colnames(pheno_tbl)
      )
    }
  }
})
