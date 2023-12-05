test_that("calculate_perf_metric() works", {
  
  set.seed(134)
  n <- 10

  dir <- withr::local_tempdir()
  n_fold <- 3
  generate_mock_data(
    n_samples = n,
    n_genes = 5,
    n_na_in_pheno = 3,
    to_csv = dir
  )

  data_spec <- DataSpec(
    name = "Mock et al. (2023)",
    directory = dir
  )
  model_spec_1 <- ModelSpec(
    name = "cox-zerosum",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold),
    response_type = "survival_censored",
    base_dir = dir
  )
  model_spec_2 <- ModelSpec(
    name = "cox-lasso",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold),
    response_type = "binary",
    base_dir = dir
  )
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(dir, "perf.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    benchmark = "ipi"
  )
  data <- read(data_spec)
  prepare_and_fit(
    expr_mat = data[["expr_mat"]],
    pheno_tbl = data[["pheno_tbl"]],
    data_spec = data_spec,
    model_spec = list(model_spec_1, model_spec_2)
  )
  expect_silent(assess_model(
    data_spec = data_spec,
    model_spec = model_spec_1,
    perf_plot_spec = perf_plot_spec
  ))
  expect_silent(assess_model(
    data_spec = data_spec,
    model_spec = model_spec_2,
    perf_plot_spec = perf_plot_spec
  ))
})
