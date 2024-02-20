test_that("plot_2d_metric() works", {

  set.seed(134)

  n_row <- 50

  dir <- withr::local_tempdir()
  ass_spec_2d <- AssSpec2d(
    file = file.path(dir, "test.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    pivot_time_cutoff = 2.,
    benchmark = "bm",
    show_plots = FALSE,
    title = "this title",
    x_lab = "this x lab",
    y_lab = "that y lab",
    xlim = c(0, .9),
    smooth_method = "loess",
    smooth_benchmark = TRUE,
    smooth_se = FALSE,
    scale_x = "identity",
    hline = list(yintercept = 0.5, linetype = "dashed"),
    text = list(ggplot2::aes(x = 0.5, y = 0.5, label = "this text")),
    colors = c("black", "blue", "yellow", "red"),
    theme = ggplot2::theme_dark() + 
      ggplot2::theme(plot.background = ggplot2::element_rect(grDevices::rgb(1, 1, 1, .98)))
  )
  ass_spec_2d$data <- tibble::tibble(
    rpp = runif(n_row),
    prec = runif(n_row),
    cutoff = sample(0:5, n_row, replace = TRUE),
    model = sample(c("model1", "model2", "bm"), n_row, replace = TRUE)
  )

  expect_no_error(
    plot_2d_metric(
      ass_spec_2d = ass_spec_2d,
      quiet = TRUE
    )
  )
  ass_spec_2d$benchmark <- NULL
  ass_spec_2d$hline <- NULL
  ass_spec_2d$text <- NULL
  ass_spec_2d$vline <- list(xintercept = 0.5, color = "blue")
  expect_no_error(
    plot_2d_metric(
      ass_spec_2d = ass_spec_2d,
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
  ass_spec_2d <- AssSpec2d(
    file = file.path(dir, "scores.pdf"),
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
      ass_spec_2d = ass_spec_2d,
      quiet = TRUE
    )
  )
  expect_true(all(
    file.exists(file.path(dir, "scores.pdf")),
    file.exists(file.path(dir, "scores.csv"))
  ))
})
