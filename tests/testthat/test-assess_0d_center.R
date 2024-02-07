test_that("assess_0d_center() works", {

  set.seed(452)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  base_dir <- withr::local_tempdir("data")
  data_dir <- file.path(base_dir, "data/mock")
  dir.create(data_dir, recursive = TRUE)
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = data_dir
  )
  data_spec <- DataSpec(
    name = "Mock et al. (2023)", 
    directory = data_dir, 
    train_prop = .7,
    benchmark_col = "ipi"
  )

  model_dir <- file.path(base_dir, "models")
  dir.create(model_dir)
  model_spec_1 <- ModelSpec(
    name = "cox-zerosum",
    directory = file.path(model_dir, "cox"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    response_type = "survival_censored",
    fit_file = "model1.rds"
  )
  model_spec_2 <- ModelSpec(
    name = "binomial-zerosum",
    directory = file.path(model_dir, "logistic"),
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    fit_file = "model2.rds"
  )
  model_spec_list <- list(model_spec_1, model_spec_2)

  training_camp(
    data_spec = data_spec,
    model_spec_list = model_spec_list,
    quiet = TRUE
  )

  res_dir <- file.path(base_dir, "results")
  ass_spec_0d <- AssSpec0d(
    metric = "get_auc",
    pivot_time_cutoff = 2,
    file = file.path(model_dir, "auc.csv")
  )

  auc_tbl <- assess_0d_center(
    model_spec_list = model_spec_list,
    data_spec = data_spec,
    ass_spec_0d = ass_spec_0d,
    cohort = c("train", "test"),
    quiet = TRUE
  )
  expect_true(
    file.exists(file.path(model_dir, "auc.csv")) &&
    file.exists(file.path(res_dir, "auc.csv"))
  )
  expect_equal(colnames(auc_tbl), c("model", "cohort", "cutoff", "mean", "sd", "min", "max"))
  expect_equal(nrow(auc_tbl), 2 * (1+2))
})
