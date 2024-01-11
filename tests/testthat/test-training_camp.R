test_that("training_camp() works", {

  set.seed(342)

  n_samples <- 20
  n_genes <- 3
  n_na_in_pheno <- 0
  n_fold <- 2
  lambda <- 1

  dir <- withr::local_tempdir()
  generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = dir
  )
  data_spec <- DataSpec(
    name = "mock",
    directory = dir,
    train_prop = .66,
    pivot_time_cutoff = 2.0
  )
  model_spec_1 <- ModelSpec(
    name = "model1",
    directory = file.path(dir, "model1"),
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 1.5,
    optional_fitter_args = list(family = "cox", nfolds = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "survival_censored"
  )
  model_spec_2 <- model_spec_1
  model_spec_2$split_index <- 1:2
  model_spec_2$time_cutoffs <- c(1.5, 2.)
  
  expect_no_error(
    training_camp(
      data_spec = data_spec,
      model_spec_list = list(model_spec_1, model_spec_2)
    )
  )
})
