test_that("plot_2d_metric() works", {

  set.seed(134)

  n_row <- 50

  dir <- withr::local_tempdir()
  ass2d <- Ass2d$new(
    x_metric = "rpp",
    y_metric = "prec",
    x_lab = "this x lab",
    y_lab = "that y lab",
    scale_x = "identity",
    colors = c("black", "blue", "yellow", "red"),
    theme = ggplot2::theme_dark() + 
      ggplot2::theme(
        plot.background = ggplot2::element_rect(grDevices::rgb(1, 1, 1, .98))
        ),
    dpi = 200
  )
  tbl <- tibble::tibble(
    rpp = runif(n_row),
    prec = runif(n_row),
    cutoff = sample(0:5, n_row, replace = TRUE),
    model = "model1"
  )

  plt <- plot_2d_metric(
    tbl = tbl,
    ass2d = ass2d,
    file = NULL,
    fellow_csv = FALSE,
    quiet = TRUE
  )
  # print(plt)

  tbl[["model"]] <- sample(c("model1", "model2"), n_row, replace = TRUE)
  plt <- plot_2d_metric(
    tbl = tbl,
    ass2d = ass2d,
    file = file.path(dir, "plt.jpeg"),
    fellow_csv = TRUE,
    quiet = TRUE
  )
  # print(plt)
  expect_true(file.exists(file.path(dir, "plt.jpeg")))
  expect_true(file.exists(file.path(dir, "plt.csv")))
})