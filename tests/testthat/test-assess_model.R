test_that("calculate_perf_metric() works", {
  
  set.seed(134)
  n <- 10

  dir <- withr::local_tempdir()
  n_fold <- 3
  data <- generate_mock_data(
    n_samples = n,
    n_genes = 5,
    n_na_in_pheno = 3
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]

  data_spec <- DataSpec(
    name = "Mock et al. (2023)"
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
  prepare_and_fit(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = list(model_spec_1, model_spec_2)
  )

  expect_silent(assess_model(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_1,
    perf_plot_spec = perf_plot_spec
  ))
  expect_silent(assess_model(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_2,
    perf_plot_spec = perf_plot_spec
  ))
})
