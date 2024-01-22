
test_that("plot_2d_metric() works", {

  set.seed(134)

  n_row <- 30

  dir <- withr::local_tempdir()
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(dir, "test.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    pivot_time_cutoff = 2.,
    benchmark = "bm",
    show_plots = TRUE,
    title = "this title",
    x_lab = "this x lab",
    y_lab = "that y lab",
    xlim = c(0, .9),
    smooth_method = "loess",
    smooth_benchmark = TRUE,
    scale_x = "log10",
    hline = list(yintercept = 0.5, linetype = "dashed"),
    text = list(ggplot2::aes(x = 0.5, y = 0.5, label = "this text"))
  )
  perf_plot_spec$data <- tibble::tibble(
    rpp = runif(n_row),
    prec = runif(n_row),
    model = sample(c("model", "bm"), n_row, replace = TRUE)
  )

  expect_no_error(
    plot_2d_metric(
      perf_plot_spec = perf_plot_spec,
      quiet = TRUE
    )
  )
  perf_plot_spec$benchmark <- NULL
  perf_plot_spec$hline <- NULL
  perf_plot_spec$text <- NULL
  perf_plot_spec$vline <- list(xintercept = 0.5, color = "blue")
  expect_no_error(
    plot_2d_metric(
      perf_plot_spec = perf_plot_spec,
      quiet = TRUE
    )
  )
})


test_that("plot_risk_scores() works", {
  
  set.seed(134)

  n_samples <- 5
  split_index <- 1:2

  l <- apb(n_samples, split_index)
  actual <- l[[1]]
  predicted <- l[[2]]
  dir <- withr::local_tempdir()
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(dir, "scores.pdf"),
    x_metric = "rank",
    y_metric = "scores",
    title = "this title",
    fellow_csv = TRUE,
    show_plots = FALSE
  )

  expect_no_error(
    plot_risk_scores(
      predicted = predicted,
      actual = actual,
      perf_plot_spec = perf_plot_spec,
      quiet = TRUE
    )
  )
  expect_true(all(
    file.exists(file.path(dir, "scores.pdf")),
    file.exists(file.path(dir, "scores.csv"))
  ))
})
