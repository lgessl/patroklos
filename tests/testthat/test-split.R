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
  data <- Data$new(
    name = "mock",
    directory = dir,
    train_prop = 0.8,
    pivot_time_cutoff = 2.0,
    time_to_event_col = "pfs_years",
    cohort = "train"
  )
  model_1 <- Model$new(
    name = "dummy1",
    directory = dir,
    fitter = ptk_zerosum,
    split_index = 1:2,
    time_cutoffs = 2,
    val_error_fun = neg_roc_auc
  )
  model_2 <- Model$new(
    name = "dummy2",
    directory = dir,
    fitter = ptk_zerosum,
    split_index = 1,
    time_cutoffs = 2,
    val_error_fun = neg_roc_auc
  )

  for(pivot_time_cutoff in list(NULL, 2)){
    data$pivot_time_cutoff <- pivot_time_cutoff
    new_pheno_tbl <- ensure_splits(
      pheno_tbl = pheno_tbl,
      data = data,
      model_list = list(model_1, model_2)
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
