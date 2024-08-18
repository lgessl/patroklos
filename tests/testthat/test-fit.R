test_that("Model$fit() works", {

  set.seed(4352)

  n_samples <- 50
  n_genes <- 5
  n_na_in_pheno <- 2
  n_fold <- 3
  lambda <- c(1, 2)

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  data$imputer <- mean_impute
  dir <- withr::local_tempdir()

  # Survival response
  model <- Model$new(
    name = "cox",
    directory = file.path(dir, "model1"),
    fitter = ptk_zerosum,
    time_cutoffs = 2,
    val_error_fun = neg_prec_with_prev_greater(0.15),
    hyperparams = list(family = "cox", nFold = n_fold, lambda = lambda, 
      zeroSum = FALSE)
  )

  model$fit(data, quiet = TRUE)
  expect_s3_class(model$fit_obj, "ptk_zerosum")
  expect_equal(nrow(model$fit_obj$coef[[1]]), n_genes+1)
  expect_true(file.exists(file.path(dir, "model1", "model.rds")))

  model$fit(data, quiet = TRUE)

  model$directory <- file.path(dir, "model2")
  model$include_from_discrete_pheno <- "discrete_var" # 3 levels
  model$include_from_continuous_pheno <- "continuous_var"
  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fit_obj$coef[[1]]), n_genes+1+2+1)
  expect_s3_class(model$fit_obj, "ptk_zerosum")

  # Binary response
  colnames(data$pheno_tbl)[1] <- "sample"
  colnames(data$pheno_tbl)[3] <- "pfs"
  data$patient_id_col <- "sample"
  data$time_to_event_col <- "pfs"
  model <- Model$new(
    name = "logistic",
    directory = file.path(dir, "model3"),
    fitter = ptk_zerosum,
    time_cutoffs = 2,
    val_error_fun = neg_roc_auc,
    hyperparams = list(family = "binomial", nFold = n_fold, lambda = lambda, 
      zeroSum = FALSE),
    include_from_discrete_pheno = c("ipi_age", "abc_gcb"),
    combine_n_max_categorical_features = 1:2,
    combined_feature_min_positive_ratio = 0
  )

  model$fit(data, quiet = TRUE)
  expect_true(nrow(model$fit_obj$coef[[1]]) %in% (1 + n_genes + c(3, 5))) 
  expect_s3_class(model$fit_obj, "ptk_zerosum")

  model$directory <- file.path(dir, "model4")
  model$include_from_discrete_pheno <- "discrete_var" # 3 levels
  model$include_from_continuous_pheno <- "continuous_var"
  model$include_expr <- FALSE
  model$fit(data, quiet = TRUE)
  expect_equal(nrow(model$fit_obj$coef[[1]]), 1+2+1)
  expect_true(file.exists(file.path(dir, "model4", "model.rds")))

  # nested_cv
  model <- Model$new(
      name = "RF pseudo OOB",
      fitter = long_nestor,
      directory = file.path(dir, "model5"),
      time_cutoffs = Inf,
      val_error_fun = error_rate,
      hyperparams = list(
        fitter1 = ptk_zerosum,
        fitter2 = multitune(ptk_ranger),
        hyperparams1 = list(
          family = "cox",
          alpha = 1,
          zeroSum = FALSE,
          lambdaSteps = 2,
          standardize = TRUE,
          nFold = n_fold
        ),
        hyperparams2 = list(
          num.trees = 100, 
          rel_mtry = FALSE,
          mtry = 2:3,
          min.node.size = 3, 
          classification = TRUE
        )
      ),
      include_from_continuous_pheno = "continuous_var",
      include_from_discrete_pheno = "discrete_var"
  )
  model$fit(data, quiet = TRUE)

  model <- Model$new(
    name = "one-column x",
    fitter = ptk_zerosum,
    directory = file.path(dir, "model6"),
    time_cutoffs = 2,
    val_error_fun = neg_prec_with_prev_greater(0.15),
    hyperparams = list(family = "binomial", nFold = n_fold, lambda = 0,
      zeroSum = FALSE),
    include_from_continuous_pheno = "continuous_var",
    include_expr = FALSE
  )
  model$fit(data, quiet = TRUE)

  model$name <- "new shell name"
  model$fit(data, quiet = TRUE, update_model_shell = TRUE)
  model <- readRDS(file.path(model$directory, model$file))
  expect_equal(model$name, "new shell name")
})