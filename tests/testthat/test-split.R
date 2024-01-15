test_that("split_dataset works", {

  set.seed(234)

  n_samples <- 100
  n_genes <- 1

  dir <- withr::local_tempdir()
  pheno_tbl <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = 0,
    split_index = NULL
  )[["pheno_tbl"]]
  pheno_tbl[["split_2"]] <- "test" # One split already there
  data_spec <- DataSpec(
    name = "mock",
    directory = dir,
    train_prop = 0.8,
    pivot_time_cutoff = 2.0
  )
  model_spec_1 <- ModelSpec(
    name = "dummy1",
    directory = dir,
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2
  )
  model_spec_2 <- ModelSpec(
    name = "dummy2",
    directory = dir,
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2
  )

  for(pivot_time_cutoff in list(NULL, 2)){
    data_spec$pivot_time_cutoff <- pivot_time_cutoff
    new_pheno_tbl <- ensure_splits(
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec_list = list(model_spec_1, model_spec_2)
    )
    expect_equal(
      nrow(new_pheno_tbl),
      nrow(pheno_tbl)
    )
    expect_equal(
      ncol(new_pheno_tbl),
      ncol(pheno_tbl) + 1
    )
    expect_true("split_1" %in% colnames(new_pheno_tbl))
  }
})
