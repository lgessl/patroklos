test_that("Model$fit() works", {

  set.seed(4352)

  n_samples <- 30
  n_genes <- 5
  n_na_in_pheno <- 2
  n_fold <- 3
  lambda <- 1

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    split_index = 1:2
  )
  dir <- withr::local_tempdir()

  # Survival response
  model <- Model$new(
    name = "cox",
    directory = file.path(dir, "model1"),
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2,
    optional_fitter_args = list(family = "cox", nfolds = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "survival_censored"
  )

  model$fit(data, quiet = TRUE)
  expect_equal(length(model$fits), 1)
  expect_equal(nrow(model$fits[[1]]$coef[[1]]), n_genes+1)
  expect_true(file.exists(file.path(dir, "model1", "models.rds")))

  model$split_index <- 1:2
  model$fit(data, quiet = TRUE)

  model$directory <- file.path(dir, "model2")
  model$split_index <- 1:2
  model$include_from_discrete_pheno <- "discrete_var" # 3 levels
  model$include_from_continuous_pheno <- "continuous_var"
  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fits[[1]]$coef[[1]]), n_genes+1+2+1)
  expect_equal(length(model$fits), length(model$split_index))

  # Binary response
  colnames(data$pheno_tbl)[1] <- "sample"
  colnames(data$pheno_tbl)[3] <- "pfs"
  data$patient_id_col <- "sample"
  data$time_to_event_col <- "pfs"
  model <- Model$new(
    name = "logistic",
    directory = file.path(dir, "model3"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2,
    optional_fitter_args = list(family = "binomial", nfolds = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "binary"
  )

  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fits[[1]]$coef[[1]]), n_genes+1)
  expect_s3_class(model$fits[[2]], "zeroSum")

  model$directory <- file.path(dir, "model4")
  model$split_index <- 1
  model$include_from_discrete_pheno <- "discrete_var" # 3 levels
  model$include_from_continuous_pheno <- "continuous_var"
  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fits[[1]]$coef[[1]]), n_genes+1+2+1)
  expect_equal(length(model$fits), length(model$split_index))
  expect_true(file.exists(file.path(dir, "model4", "models.rds")))
})