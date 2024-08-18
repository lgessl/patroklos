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
  data$imputer <- mean_impute
  
  model_1 <- Model$new(
    name = "model1",
    directory = file.path(dir, "model1"),
    fitter = ptk_zerosum,
    time_cutoffs = c(1.5, 2),
    val_error_fun = neg_prec_with_prev_greater(0.15),
    hyperparams = list(family = "cox", nFold = n_fold, lambda = lambda, 
      zeroSum = FALSE)
  )

  model_2 <- Model$new(
      name = "RF pseudo OOB",
      fitter = long_nestor,
      directory = file.path(dir, "model2"),
      time_cutoffs = 2.,
      val_error_fun = error_rate,
      hyperparams = list(
          fitter1 = ptk_zerosum,
          fitter2 = multitune(ptk_ranger),
          hyperparams1 = list(
              family = "binomial",
              alpha = 1,
              zeroSum = FALSE,
              lambda = c(0.5, 1)
          ),
          hyperparams2 = list(
              num.trees = 15, # If chosen too small, OOB predictions can be NaN
              rel_mtry = FALSE,
              mtry = 2,
              min.node.size = 3, 
              classification = TRUE
          )
      ),
      include_from_continuous_pheno = "continuous_var",
      include_from_discrete_pheno = "discrete_var"
  )
  
  model_3 <- model_1$clone()
  model_3$directory <- file.path(dir, "model3")
  model_3$time_cutoffs <- 2
  model_3$hyperparams[["family"]] <- "binomial"

  model4 <- Model$new(
    name = "cnn projection",
    directory = file.path(dir, "cnn"),
    fitter = projection_on_feature,
    time_cutoffs = 2,
    val_error_fun = neg_prec_with_prev_greater(0.17),
    hyperparams = list(feature = "continuous_var"),
    include_from_continuous_pheno = "continuous_var",
    include_expr = FALSE
  )
  
  training_camp(
    model_list = list(model_1, model_2, model_3, model4),
    data = data,
    quiet = TRUE,
    skip_on_error = FALSE
  )
  
  expect_true(file.exists(file.path(model_1$directory, "model.rds")))
  expect_true(file.exists(file.path(model_2$directory, "model.rds")))
  expect_true(file.exists(file.path(model_3$directory, "model.rds")))
  expect_true(file.exists(file.path(model4$directory, "model.rds")))
  # expect_false(file.exists(file.path(model_4$directory, "model.rds")))
})
