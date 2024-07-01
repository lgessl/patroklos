test_that("plot_2d_metric() works", {

  set.seed(134)

  n_row <- 50

  dir <- withr::local_tempdir()
  ass2d <- Ass2d$new(
    file = file.path(dir, "test.jpeg"),
    x_metric = "rpp",
    y_metric = "prec",
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
      ggplot2::theme(
        plot.background = ggplot2::element_rect(grDevices::rgb(1, 1, 1, .98)),
        text = ggplot2::element_text(family = "Courier")
        ),
    dpi = 200
  )
  ass2d$data <- tibble::tibble(
    rpp = runif(n_row),
    prec = runif(n_row),
    cutoff = sample(0:5, n_row, replace = TRUE),
    model = sample(c("model1", "model2", "bm"), n_row, replace = TRUE)
  )

  expect_no_error(
    plot_2d_metric(
      ass2d = ass2d,
      quiet = TRUE
    )
  )
  ass2d$benchmark <- NULL
  ass2d$hline <- NULL
  ass2d$text <- NULL
  ass2d$vline <- list(xintercept = 0.5, color = "blue")
  expect_no_error(
    plot_2d_metric(
      ass2d = ass2d,
      quiet = TRUE
    )
  )
})