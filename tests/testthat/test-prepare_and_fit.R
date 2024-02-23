test_that("prepare_and_fit() works", {

  set.seed(4352)

  n_samples <- 20
  n_genes <- 5
  n_na_in_pheno <- 0
  n_fold <- 3
  lambda <- 1
  split_index <- 1:2

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    split_index = split_index
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]
  pheno_tbl[["split_1"]] <- sample(c("train", "test"), n_samples, replace = TRUE)
  pheno_tbl[["split_2"]] <- sample(c("train", "test"), n_samples, replace = TRUE)
  data <- DataSpec(
    name = "Mock et al. (2023)", 
    directory = "mock", 
    train_prop = 0.8,
    cohort = "train"
  )
  dir <- withr::local_tempdir()

  model <- ModelSpec(
    name = "cox-vanilla",
    directory = file.path(dir, "model1"),
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2,
    optional_fitter_args = list(family = "cox", nfolds = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "survival_censored",
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL
  )

  fits <- prepare_and_fit(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data = data,
    model = model,
    quiet = TRUE
  )
  expect_equal(length(fits), 2)

  model$split_index <- split_index
  fits <- prepare_and_fit(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data = data,
    model = model,
    quiet = TRUE
  )
  expect_equal(length(fits), length(split_index)+1)
})