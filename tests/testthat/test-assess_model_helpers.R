test_that("plot_perf_metric() works", {

  set.seed(134)

  n_samples <- 5

  perf_tbl <- tibble::tibble(
    rpp = runif(2*n_samples),
    prec = runif(2*n_samples),
    cutoff = runif(2*n_samples),
    model = sample(c("model1", "model2"), 2*n_samples, replace = TRUE)
  )
  dir <- withr::local_tempdir()
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(dir, "test.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    show_plots = FALSE,
    title = "this title",
    x_lab = "this x lab",
    y_lab = "that y lab"
  )
  perf_plot_spec$data <- perf_tbl

  expect_silent(
    plot_perf_metric(
      perf_plot_spec = perf_plot_spec,
      quiet = TRUE
    )
  )
})

test_that("calculate_perf_metric() works", {
    
    set.seed(134)

    n_samples <- 5

    predicted <- rnorm(n_samples)
    actual <- sample(c(0, 1), n_samples, replace = TRUE)
    perf_plot_spec <- PerfPlotSpec(
      fname = "test.pdf",
      x_metric = "rpp",
      y_metric = "prec"
    )
  
    expect_silent(
      perf_plot_spec <- calculate_perf_metric(
        predicted = predicted,
        actual = actual,
        perf_plot_spec = perf_plot_spec
      )
    )
    perf_tbl <- perf_plot_spec$data
    expect_equal(names(perf_tbl), c("rpp", "prec", "cutoff"))
    expect_s3_class(perf_tbl, "tbl_df")
    expect_equal(dim(perf_tbl), c(n_samples, 3))
})

test_that("add_benchmark_perf_metric() works", {

  set.seed(134)

  n_samples <- 20

  pheno_tbl <- generate_mock_data(
    n_samples = n_samples,
    n_genes = 1,
    n_na_in_pheno = 0
  )[["pheno_tbl"]]
  pheno_tbl[["ipi"]][1] <- NA
  data_spec <- DataSpec(name = "Mock et al. (2023)")
  model_spec <- ModelSpec(
    name = "cox-zerosum",
    fitter = zeroSum::zeroSum
  )
  perf_plot_spec <- PerfPlotSpec(
    fname = "test.pdf",
    x_metric = "rpp",
    y_metric = "prec"
  )

  perf_plot_spec <- add_benchmark_perf_metric(
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    perf_plot_spec = perf_plot_spec,
    model_spec = model_spec
  )
  perf_tbl <- perf_plot_spec$bm_data

  expect_equal(names(perf_tbl), c("rpp", "prec", "cutoff"))
  expect_s3_class(perf_tbl, "tbl_df")
  expect_lt(nrow(perf_tbl), n_samples)
  expect_equal(ncol(perf_tbl), 3)  
})

test_that("plot_scores() works", {
  
    set.seed(134)
  
    n_samples <- 5
  
    predicted <- rnorm(n_samples)
    actual <- sample(c(0, 1), n_samples, replace = TRUE)
    dir <- withr::local_tempdir()
    perf_plot_spec <- PerfPlotSpec(
      fname = file.path(dir, "scores.pdf"),
      x_metric = "rank",
      y_metric = "scores",
      title = "this title",
      show_plots = FALSE
    )

    expect_silent(
      plot_scores(
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
