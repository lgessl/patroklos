test_that("AssScalar$new() works", {

    ass_scalar <- AssScalar$new(
        metrics = c("auc", "accuracy", "prevalence", "precision_ci_ll", "n_true"),
        prev_range = c(0.15, 0.30),
        confidence_level = 0.9,
        benchmark = "ipi",
        file = "some file",
        round_digits = 4
    )
    expect_equal(ass_scalar$metrics, c("auc", "accuracy", "prevalence", "precision_ci_ll", "n_true"))
    expect_equal(ass_scalar$prev_range, c(0.15, 0.30))
    expect_equal(ass_scalar$confidence_level, 0.9)
    expect_equal(ass_scalar$benchmark, "ipi")
    expect_equal(ass_scalar$round_digits, 4)
    expect_equal(ass_scalar$file, "some file")

    expect_error(AssScalar$new(metrics = c("auc", "blabla")))
})

test_that("AssScalar$assess() works", {

    set.seed(342)
    n_samples <- 50
    n_genes <- 5
    n_fold <- 1
    metrics <- c("auc", "accuracy", "precision", "logrank")
    
    model_dir <- withr::local_tempdir()
    data <- generate_mock_data(
        n_samples = n_samples,
        n_genes = n_genes,
        n_na_in_pheno = 1
    )
    model <- Model$new(
      name = "rf",
      directory = file.path(model_dir, "rf"),
      fitter = hypertune(ptk_ranger, select = TRUE),
      hyperparams = list(rel_mtry = TRUE, mtry = 0.98, num.trees = 100, 
        min.node.size = 3, classification = TRUE),
      time_cutoffs = 2,
      val_error_fun = error_rate,
      file = "model.rds"
    )
    data$cohort <- "train"
    model$fit(data, quiet = TRUE)
    ass_scalar <- AssScalar$new(
      metrics = metrics,
      benchmark = "ipi",
      round_digits = 4
    )

    data$cohort <- "test"
    res <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(is.numeric(res))
    expect_true(is.vector(res))
    expect_equal(length(res), length(metrics))

    ass_scalar$metrics <- c("precision", "prevalence")
    ass_scalar$prev_range <- c(0.15, 1)
    res <- ass_scalar$assess(data, model, quiet = TRUE)

    ass_scalar$prev_range <- c(0.15, 1)
    model <- Model$new(
      name = "log",
      directory = file.path(model_dir, "log"),
      fitter = ptk_zerosum,
      hyperparams = list(family = "gaussian", lambdaSteps = 1, zeroSum = FALSE, 
        nFold = 1),
      time_cutoffs = 2,
      val_error_fun = neg_prec_with_prev_greater(0.15),
    )
    model$fit(data, quiet = TRUE)
    res <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(res[2] >= 0.15 && res[2] <= 1)

    ass_scalar$prev_range <- c(0.1000, 0.10000)
    res <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(all(is.na(res)))
})

test_that("AssScalar$assess_center() works", {

  set.seed(452)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 2
  multiple_metrics <- c("precision", "accuracy", "n_true", "perc_true", "n_samples", 
    "logrank", "threshold", "prevalence", "precision_ci_ll", "hr", "hr_ci_ll", "hr_ci_ul", "hr_p")
  single_metric <- "auc"

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = file.path(dir, "data")
  )
  model1 <- Model$new(
    name = "cox",
    directory = file.path(dir, "models/cox"),
    fitter = ptk_zerosum,
    hyperparams = list(family = "cox", nFold = n_fold, 
      lambdaSteps = 2, zeroSum = FALSE),
    time_cutoffs = c(1.5, 2),
    val_error_fun = neg_prec_with_prev_greater(0.17),
    file = "model2.rds"
  )
  model2 <- Model$new(
    name = "rf",
    directory = file.path(dir, "models/rf"),
    fitter = hypertune(ptk_ranger),
    hyperparams = list(rel_mtry = FALSE, mtry = 2, num.trees = 100, min.node.size = 3,
      classification = TRUE),
    time_cutoffs = 2,
    val_error_fun = error_rate,
    file = "model3.rds",
    include_from_discrete_pheno = c("discrete_var", "abc_gcb", "ipi_age"),
    combine_n_max_categorical_features = 3,
    combined_feature_min_positive_ratio = 0.04
  )
  model3 <- Model$new(
    name = "greedy nestor",
    directory = file.path(dir, "models/nestor"),
    fitter = greedy_nestor,
    time_cutoffs = 2,
    hyperparams = list(
      model1 = model1,
      fitter2 = ptk_zerosum,
      hyperparams2 = list(family = "gaussian", nFold = 2, lambdaSteps = 2, zeroSum = FALSE)
    ),
    include_from_discrete_pheno = c("discrete_var", "abc_gcb"),
    val_error_fun = neg_prec_with_prev_greater(0.20),
    combine_n_max_categorical_features = 2,
    combined_feature_min_positive_ratio = 0.05
  )
  model4 <- Model$new(
    name = "greedy nestor rf",
    directory = file.path(dir, "models/nestor_rf"),
    fitter = greedy_nestor,
    time_cutoffs = 2,
    hyperparams = list(
      model1 = model1,
      fitter2 = hypertune(ptk_ranger),
      hyperparams2 = list(rel_mtry = FALSE, num.trees = 100, min.node.size = c(4, 5), 
        classification = TRUE)
    ),
    include_from_discrete_pheno = c("discrete_var", "abc_gcb"),
    val_error_fun = error_rate,
    combine_n_max_categorical_features = 2,
    combined_feature_min_positive_ratio = 0.05
  )
  model5 <- Model$new(
    name = "nemo",
    directory = file.path(dir, "nemo"),
    fitter = ptk_zerosum,
    time_cutoffs = 2,
    val_error_fun = neg_prec_with_prev_greater(0.15),
    hyperparams = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambdaSteps = 1, zeroSum = FALSE)
  )

  model1_list <- list(model1, model3)
  model2_list <- list(model1, model2, model4)
  model_list <- list(model1, model2, model3, model4)
  training_camp(model_list, data, quiet = TRUE, skip_on_error = FALSE) 
  model_list <- c(model_list, list(model5))
  ass_scalar <- AssScalar$new(
    metrics = "auc", # just for the moment
    file = file.path(dir, "models/eval.csv"),
    prev_range = c(0.15, 1)
  )
  dir.create(dir, "results")
  
  data$cohort <- "val_predict"
  ass_scalar$metrics <- multiple_metrics
  eval_tbl <- ass_scalar$assess_center(data, model2_list, quiet = TRUE)
  expect_true(file.exists(file.path(dir, "models/eval.csv")))
  expect_equal(nrow(eval_tbl), length(model2_list))  
  expect_equal(colnames(eval_tbl), c("model", "cohort", multiple_metrics))
  expect_true(is.numeric(as.matrix(eval_tbl[, 4:ncol(eval_tbl)])))
})

test_that("prepend_to_filename works", {

  ass_scalar1 <- AssScalar$new(metrics = "accuracy", file = "file1")
  ass_scalar2 <- ass_scalar1$clone()
  ass_scalar2$file <- "file2"
  as_list <- list(ass_scalar1, ass_scalar2)
  # Prepend prefix to the file attributes
  prepend_to_filename(as_list, "prefix")

  # Check that the file attributes have been correctly prepended
  expect_equal(as_list[[1]]$file, "prefix/file1")
  expect_equal(ass_scalar1$file, "prefix/file1")
  expect_equal(ass_scalar2$file, "prefix/file2")
})