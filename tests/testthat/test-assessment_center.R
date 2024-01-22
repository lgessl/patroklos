test_that("assessment_center() works", {

  set.seed(4352)

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
    train_prop = .7
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
    fit_fname = "model1.rds"
  )
  model_spec_2 <- ModelSpec(
    name = "binomial-zerosum",
    directory = file.path(model_dir, "logistic"),
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    fit_fname = "model2.rds"
  )
  model_spec_list <- list(model_spec_1, model_spec_2)

  training_camp(
    data_spec = data_spec,
    model_spec_list = model_spec_list
  )

  res_dir <- file.path(base_dir, "results")
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(model_dir, "perf_plot.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    show_plots = FALSE,
    smooth_method = "loess",
    pivot_time_cutoff = 2.,
    fellow_csv = TRUE,
    scores_plot = TRUE
  )

  assessment_center(
      model_spec_list = model_spec_list,
      data_spec = data_spec,
      perf_plot_spec = perf_plot_spec
  )
  expect_true(file.exists(perf_plot_spec$fname))
  expect_true(file.exists(file.path(model_dir, "logistic/1-5/rpp_vs_prec.pdf")))
  expect_true(file.exists(file.path(res_dir, "logistic/2/scores.pdf")))
  expect_true(file.exists(file.path(res_dir, "cox/2/rpp_vs_prec.csv")))

  perf_plot_spec$y_metric <- "logrank"
  perf_plot_spec$fname <- file.path(model_dir, "logrank.pdf")
  perf_plot_spec$benchmark <- NULL
  perf_plot_spec$scores_plot <- FALSE
  assessment_center(
    model_spec_list = model_spec_list,
    data_spec = data_spec,
    perf_plot_spec = perf_plot_spec,
    comparison_plot = FALSE
  )
})
