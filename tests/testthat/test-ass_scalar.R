test_that("AssScalar$new() works", {

    ass_scalar <- AssScalar$new(
        metrics = c("auc", "accuracy", "prevalence", "n_true"),
        prev_range = c(0.15, 0.30),
        benchmark = "ipi",
        file = "some file",
        round_digits = 4
    )
    expect_equal(ass_scalar$metrics, c("auc", "accuracy", "prevalence", "n_true"))
    expect_equal(ass_scalar$prev_range, c(0.15, 0.30))
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
    split_index <- 1:2
    metrics <- c("auc", "accuracy", "precision", "logrank")
    
    model_dir <- withr::local_tempdir()
    data <- generate_mock_data(
        n_samples = n_samples,
        n_genes = n_genes,
        n_na_in_pheno = 1,
        split_index = split_index
    )
    model <- Model$new(
      name = "rf",
      directory = file.path(model_dir, "rf"),
      fitter = hypertune(ptk_ranger, select = TRUE),
      hyperparams = list(rel_mtry = TRUE, mtry = 0.98, num.trees = 100, 
        min.node.size = 3, classification = TRUE),
      split_index = 1:2,
      time_cutoffs = 2,
      val_error_fun = error_rate,
      fit_file = "model.rds"
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
    expect_true(is.matrix(res))
    expect_equal(dim(res), c(length(split_index), length(metrics)))

    ass_scalar$metrics <- c("precision", "prevalence")
    ass_scalar$prev_range <- c(0.15, 0.30)
    res <- ass_scalar$assess(data, model, quiet = TRUE)

    ass_scalar$prev_range <- c(0.15, 0.30)
    model <- Model$new(
      name = "log",
      directory = file.path(model_dir, "log"),
      fitter = ptk_zerosum,
      hyperparams = list(family = "gaussian", lambda = 0, zeroSum = FALSE, 
        nFold = 1),
      split_index = 1:2,
      time_cutoffs = 2,
      val_error_fun = neg_prec_with_prev_greater(0.15),
    )
    model$fit(data, quiet = TRUE)
    res <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(all(res[, 2] >= 0.15 & res[, 2] <= 0.30))

    ass_scalar$prev_range <- c(0.1000, 0.10000)
    res <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(all(is.na(res)))
})

test_that("AssScalar$assess_center() works", {

  set.seed(452)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1
  multiple_metrics <- c("precision", "accuracy", "n_true", "perc_true", "n_samples", 
    "logrank", "threshold")
  single_metric <- "auc"

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = file.path(dir, "data")
  )
  model1 <- Model$new(
    name = "binomial-zerosum",
    directory = file.path(dir, "models/logistic"),
    fitter = ptk_zerosum,
    hyperparams = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
    val_error_fun = neg_binomial_log_likelihood,
    include_from_discrete_pheno = c("discrete_var", "abc_gcb", "ipi_age"),
    fit_file = "model2.rds",
    combine_n_max_categorical_features = 3,
    combined_feature_min_positive_ratio = 0.04
  )
  model2 <- Model$new(
    name = "rf",
    directory = file.path(dir, "models/rf"),
    fitter = hypertune(ptk_ranger),
    hyperparams = list(rel_mtry = FALSE, mtry = 2, num.trees = 100, min.node.size = 3,
      classification = TRUE),
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
    val_error_fun = error_rate,
    fit_file = "model3.rds",
    include_from_discrete_pheno = "discrete_var"
  )
  model1_list <- list(model1, model1)
  model2_list <- list(model1, model2)
  model_list <- list(model1, model2)
  training_camp(model_list, data, quiet = TRUE, skip_on_error = FALSE) 
  ass_scalar <- AssScalar$new(
    metrics = "auc", # just for the moment
    file = file.path(dir, "models/eval.csv"),
    prev_range = c(0.15, 0.45)
  )
  dir.create(dir, "results")
  
  # Case 1: one metric
  data$cohort <- "test"
  ass_scalar$metrics <- single_metric 
  eval_tbl <- ass_scalar$assess_center(data, model1_list, quiet = TRUE)
  expect_equal(nrow(eval_tbl), length(model1_list))  
  expect_equal(colnames(eval_tbl), c("model", "cohort", "mean", "sd", "min", "max"))
  csv_path <- ass_scalar$file
  first_2 <- readr::read_lines(csv_path, n_max = 2)
  expect_equal(first_2[1], "# mock | pfs_years < 2")
  expect_equal(first_2[2], paste0(colnames(eval_tbl), collapse = ","))
  eval_tbl_read <- readr::read_csv(csv_path, comment = "#", show_col_types = FALSE)
  expect_equal(dim(eval_tbl), dim(eval_tbl_read))

  # Case 2: multiple metrics
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