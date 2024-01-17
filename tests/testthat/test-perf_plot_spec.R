test_that("infer_pps() works", {

  mirror <- c("models", "results")
  perf_plot_spec <- PerfPlotSpec(
    fname = "perf_plot_spec.png",
    x_metric = "x",
    y_metric = "y",
    pivot_time_cutoff = 2.,
    benchmark = "benchmark",
    fellow_csv = TRUE
  )
  perf_plot_spec$model_tree_mirror <- mirror
  data_spec <- DataSpec(
    name = "mock",
    directory = "data/mock",
    train_prop = .7,
    cohort = "train"
  )
  model_spec <- ModelSpec(
    name = "cox",
    directory = "models/cox",
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.
  )

  auto_pps <- infer_pps(
    perf_plot_spec = perf_plot_spec,
    model_spec = model_spec,
    data_spec = data_spec
  )
  expect_equal(
    auto_pps$fname,
    "models/cox/x_vs_y.png"
  )
  expect_equal(
    auto_pps$title,
    "mock train, pfs_years < 2"
  )

  data_spec$cohort <- "test"
  auto_pps <- infer_pps(
    perf_plot_spec = perf_plot_spec,
    model_spec = model_spec,
    data_spec = data_spec
  )
  expect_equal(
    auto_pps$fname,
    "results/cox/x_vs_y.png"
  )
  expect_equal(
    auto_pps$title,
    "mock test, pfs_years < 2"
  )
})
