test_that("calculate_perf_metric() works", {
  
  set.seed(134)
  n <- 10
  predicted <- rnorm(n)
  actual <- sample(c(0, 1), n, replace = TRUE)
  perf_plot_spec <- PerfPlotSpec(
    fname = "test.pdf",
    x_metric = "rpp",
    y_metric = "prec"
  )

  expect_silent(
    perf_tbl <- calculate_perf_metric(
      predicted = predicted,
      actual = actual,
      perf_plot_spec = perf_plot_spec
    )
  )
  expect_equal(names(perf_tbl), c("rpp", "prec", "cutoff"))
  expect_s3_class(perf_tbl, "tbl_df")
  testthat::expect_true(nrow(perf_tbl) <= n)
})
