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
    
    model_dir <- withr::local_tempdir()
    data <- generate_mock_data(
        n_samples = n_samples,
        n_genes = n_genes,
        n_na_in_pheno = 1,
        split_index = split_index
    )
    model <- Model$new(
        name = "logistic",
        directory = model_dir,
        fitter = zeroSum::zeroSum,
        split_index = split_index,
        time_cutoffs = 2.,
        optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
            lambda = lambda, zeroSum = FALSE),
        response_type = "binary"
    )
    ass_scalar <- AssScalar$new(
        metric = "get_auc",
        pivot_time_cutoff = 2,
        benchmark = "ipi",
        round_digits = 4
    )
    model$fit(data, quiet = TRUE)

    auc <- ass_scalar$assess(data, model, quiet = TRUE)
    expect_true(is.numeric(auc))
    expect_equal(length(auc), length(split_index))
})

test_that("AssScalar$assess_center() works", {

  set.seed(452)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = file.path(dir, "data")
  )
  model_1 <- Model$new(
    name = "cox-zerosum",
    directory = file.path(dir, "models/cox"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    response_type = "survival_censored",
    include_from_continuous_pheno = "continuous_var",
    fit_file = "model1.rds"
  )
  model_2 <- Model$new(
    name = "binomial-zerosum",
    directory = file.path(dir, "models/logistic"),
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    include_from_discrete_pheno = "discrete_var",
    fit_file = "model2.rds"
  )
  model_list <- list(model_1, model_2)
  training_camp(model_list, data, quiet = TRUE) 
  ass_scalar <- AssScalar$new(
    metric = "get_auc",
    pivot_time_cutoff = 2,
    file = file.path(dir, "models/auc.csv")
  )
  dir.create(dir, "results")

  data$cohort <- c("train", "test")
  auc_tbl <- ass_scalar$assess_center(data, model_list, quiet = TRUE)

  expect_true(
    file.exists(file.path(dir, "models/auc.csv")) &&
    file.exists(file.path(dir, "results/auc.csv"))
  )
  expect_equal(colnames(auc_tbl), c("model", "cohort", "cutoff", "mean", "sd", "min", "max"))
  expect_equal(nrow(auc_tbl), 2 * (1+2))  
})
