test_that("calculate_2d_metric() works", {
    
  set.seed(134)

  n_samples <- 20
  split_index <- 1:2

  l <- apb(n_samples, split_index)
  actual <- l[[1]]
  predicted <- l[[2]]
  benchmark <- l[[3]]
  model_spec <- ModelSpec(
    name = "mock",
    directory = "mock",
    fitter = zeroSum::zeroSum,
    split_index = split_index,
    time_cutoffs = 2.
  )
  perf_plot_spec <- PerfPlotSpec(
    fname = "test.pdf",
    x_metric = "rpp",
    y_metric = "prec"
  )

  expect_silent(
    perf_plot_spec <- calculate_2d_metric(
      actual = actual,
      predicted = predicted,
      benchmark = benchmark,
      perf_plot_spec = perf_plot_spec,
      model_spec = model_spec
    )
  )
  perf_tbl <- perf_plot_spec$data
  expect_equal(names(perf_tbl), c("rpp", "prec", "cutoff", "split", "model"))
  expect_s3_class(perf_tbl, "tbl_df")
  expect_true(all(perf_tbl[["model"]] %in% c(model_spec$name, perf_plot_spec$benchmark)))
  expect_true(perf_tbl[, 1:4] |> as.matrix() |> is.numeric() |> all())

  perf_plot_spec$benchmark <- NULL
  expect_silent(
    perf_plot_spec <- calculate_2d_metric(
      actual = actual,
      predicted = predicted,
      benchmark = NULL,
      perf_plot_spec = perf_plot_spec,
      model_spec = model_spec
    )
  )
})