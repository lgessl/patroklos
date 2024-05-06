test_that("AssScalar$new() works", {

    ass_scalar <- AssScalar$new(
        metric = "get_auc",
        pivot_time_cutoff = 2,
        lambda = 1.5,
        benchmark = "ipi",
        file = "some file",
        round_digits = 4
    )
    expect_equal(ass_scalar$metric, "get_auc")
    expect_equal(ass_scalar$pivot_time_cutoff, 2)
    expect_equal(ass_scalar$lambda, 1.5)
    expect_equal(ass_scalar$benchmark, "ipi")
    expect_equal(ass_scalar$round_digits, 4)
    expect_equal(ass_scalar$file, "some file")
})

test_that("AssScalar$assess() works", {

    set.seed(342)
    n_samples <- 50
    n_genes <- 5
    n_fold <- 1
    lambda <- 1
    split_index <- 1:2
    metric <- c("auc", "accuracy", "precision", "logrank")
    
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
      fitter = ptk_ranger,
      hyperparams = list(mtry = 2, num.trees = 100, min.node.size = 3),
      split_index = 1:2,
      time_cutoffs = 2,
      response_type = "binary",
      fit_file = "model.rds"
    )
    ass_scalar <- AssScalar$new(
      metric = metric,
      pivot_time_cutoff = 2,
      benchmark = "ipi",
      round_digits = 4
    )
    model$fit(data, quiet = TRUE)

    res <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(is.numeric(res))
    expect_true(is.matrix(res))
    expect_equal(dim(res), c(length(split_index), length(metric)))
})

test_that("AssScalar$assess_center() works", {

  set.seed(452)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1
  multiple_metric <- c("accuracy", "precision")
  single_metric <- "accuracy"

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
    fitter = zeroSum::zeroSum,
    hyperparams = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    include_from_discrete_pheno = "discrete_var",
    fit_file = "model2.rds"
  )
  model2 <- Model$new(
    name = "rf",
    directory = file.path(dir, "models/rf"),
    fitter = ptk_ranger,
    hyperparams = list(mtry = 2, num.trees = 100, min.node.size = 3),
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    fit_file = "model3.rds",
    include_from_discrete_pheno = "discrete_var"
  )
  model1_list <- list(model1, model1)
  model2_list <- list(model2, model2)
  model_list <- list(model1, model2)
  training_camp(model_list, data, quiet = TRUE) 
  ass_scalar <- AssScalar$new(
    metric = "auc", # just for the moment
    pivot_time_cutoff = 2,
    file = file.path(dir, "models/eval.csv")
  )
  dir.create(dir, "results")
  
  # Case 1: one metric
  data$cohort <- "test"
  ass_scalar$metric <- single_metric 
  eval_tbl <- ass_scalar$assess_center(data, model1_list, quiet = TRUE)
  expect_equal(nrow(eval_tbl), length(data$cohort) * (2+2))  
  expect_equal(colnames(eval_tbl), c("model", "cohort", "cutoff", "mean", 
    "sd", "min", "max"))

  # Case 2: multiple metrics
  data$cohort <- c("train", "test")
  ass_scalar$metric <- multiple_metric
  eval_tbl <- ass_scalar$assess_center(data, model2_list, quiet = TRUE)
  expect_true(
    file.exists(file.path(dir, "models/eval.csv")) &&
    file.exists(file.path(dir, "results/eval.csv"))
  )
  expect_equal(nrow(eval_tbl), length(data$cohort) * (2+2))  
  expect_equal(colnames(eval_tbl), c("model", "cohort", "cutoff", multiple_metric))
})
