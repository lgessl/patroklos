test_that("training_camp() works", {

  set.seed(342)

  n_samples <- 100 
  n_genes <- 3
  n_na_in_pheno <- 0
  n_fold <- 2
  lambda <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = dir
  )
  
  model_1 <- Model$new(
    name = "model1",
    directory = file.path(dir, "model1"),
    fitter = ptk_zerosum,
    split_index = 1,
    time_cutoffs = c(1.5, 2),
    hyperparams = list(family = "cox", nFold = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "survival_censored"
  )

  model_2 <- Model$new(
      name = "RF pseudo OOB",
      fitter = nested_pseudo_cv,
      directory = file.path(dir, "model2"),
      split_index = 1,
      time_cutoffs = 2.,
      hyperparams = list(
          fitter1 = ptk_zerosum,
          fitter2 = ptk_ranger,
          hyperparams1 = list(
              family = "binomial",
              alpha = 1,
              zeroSum = FALSE
          ),
          hyperparams2 = list(
              num.trees = 10, 
              mtry = 2,
              min.node.size = 3, 
              classification = TRUE
          ),
          oob = c(FALSE, TRUE)
      ),
      response_type = "binary",
      include_from_continuous_pheno = "continuous_var",
      include_from_discrete_pheno = "discrete_var"
  )
  model_2$split_index <- 1:2
  model_2$time_cutoffs <- 1.75 
  
  model_3 <- model_1$clone()
  model_3$directory <- file.path(dir, "model3")
  model_3$time_cutoffs <- 2
  model_3$hyperparams[["family"]] <- "binomial"
  model_3$response_type <- "binary"
  
  training_camp(
    model_list = list(model_1, model_2, model_3),
    data = data,
    quiet = TRUE
  )
  expect_true(file.exists(file.path(model_1$directory, "1-5", "models.rds")))
  expect_true(file.exists(file.path(model_2$directory, "1-75", "models.rds")))
  expect_true(file.exists(file.path(model_3$directory, "2", "models.rds")))
})
