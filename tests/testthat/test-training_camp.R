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
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = c(1.5, 2),
    hyperparams = list(family = "cox", nfolds = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    response_type = "survival_censored"
  )
  model_2 <- Model$new(
      name = "RF pseudo OOB",
      fitter = nested_cv_oob,
      directory = file.path(dir, "model2"),
      split_index = 1,
      time_cutoffs = 2.,
      hyperparams = list(
          fitter1 = zeroSum::zeroSum,
          fitter2 = ranger::ranger,
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
          n_folds = 2,
          pseudo_cv = TRUE
      ),
      response_type = "binary",
      include_from_continuous_pheno = "continuous_var",
      include_from_discrete_pheno = "discrete_var"   
  )
  # model_2 <- model_1$clone()
  model_2$split_index <- 1:2
  model_2$time_cutoffs <- 1.75 
  
  training_camp(
    model_list = list(model_1, model_2),
    data = data,
    quiet = TRUE
  )
  expect_true(file.exists(file.path(model_1$directory, "1-5", "models.rds")))
  expect_true(file.exists(file.path(model_2$directory, "1-75", "models.rds")))
})
