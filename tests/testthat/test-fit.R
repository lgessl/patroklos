test_that("Model$fit() works", {

  set.seed(4352)

  n_samples <- 30
  n_genes <- 5
  n_na_in_pheno <- 2
  n_fold <- 3
  lambda <- c(1, 2)

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
    fitter = ptk_zerosum,
    split_index = 1,
    time_cutoffs = 2,
    hyperparams = list(family = "cox", nFold = n_fold, lambda = lambda, 
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
    fitter = ptk_zerosum,
    split_index = 1:2,
    time_cutoffs = 2,
    hyperparams = list(family = "binomial", nFold = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "binary"
  )

  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fits[[1]]$coef[[1]]), n_genes+1)
  expect_s3_class(model$fits[[2]], "ptk_zerosum")

  model$directory <- file.path(dir, "model4")
  model$split_index <- 1
  model$include_from_discrete_pheno <- "discrete_var" # 3 levels
  model$include_from_continuous_pheno <- "continuous_var"
  model$include_expr <- FALSE
  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fits[[1]]$coef[[1]]), 1+2+1)
  expect_equal(length(model$fits), length(model$split_index))
  expect_true(file.exists(file.path(dir, "model4", "models.rds")))

  # nested_cv
  model <- Model$new(
      name = "RF pseudo OOB",
      fitter = nested_pseudo_cv,
      directory = file.path(dir, "model5"),
      split_index = 1,
      time_cutoffs = 2.,
      hyperparams = list(
        fitter1 = ptk_zerosum,
        fitter2 = hypertune(ptk_ranger, metric = "accuracy"),
        hyperparams1 = list(
          family = "binomial",
          alpha = 1,
          zeroSum = FALSE,
          lambda = lambda
        ),
        hyperparams2 = list(
          num.trees = 100, 
          mtry = 2:3,
          min.node.size = 3, 
          classification = TRUE
        )
      ),
      response_type = "binary",
      include_from_continuous_pheno = "continuous_var",
      include_from_discrete_pheno = "discrete_var"
  )
  model$fit(data, quiet = TRUE)
})