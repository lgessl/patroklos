test_that("Ass2d$new() works", {

  ass2d <- Ass2d$new(
    x_metric = "rpp",
    y_metric = "prec",
    pivot_time_cutoff = 2,
    lambda = 1.6,
    benchmark = "ipi",
    file = "some-file",
    ci_level = .66,
    fellow_csv = TRUE,
    scores_plot = TRUE,
    show_plots = FALSE,
    title = "title",
    x_lab = "x_lab",
    y_lab = "y_lab",
    xlim = c(0, .5),
    ylim = c(0, .6),
    smooth_method = "loess",
    smooth_benchmark = TRUE,
    smooth_se = FALSE,
    scale_x = "log10",
    scale_y = "log10",
    dpi = 250,
    units = "cm",
    text_size = 2.2,
    alpha = 1,
    width = 5,
    height = 5,
    colors = c("red", "blue")
  )
  expect_s3_class(ass2d, "Ass2d")
  expect_equal(ass2d$x_metric, "rpp")
  expect_equal(ass2d$y_metric, "prec")
  expect_equal(ass2d$pivot_time_cutoff, 2)
  expect_equal(ass2d$lambda, 1.6)
  expect_equal(ass2d$benchmark, "ipi")
  expect_equal(ass2d$file, "some-file")
  expect_equal(ass2d$ci_level, .66)
  expect_equal(ass2d$fellow_csv, TRUE)
  expect_equal(ass2d$scores_plot, TRUE)
  expect_equal(ass2d$show_plots, FALSE)
  expect_equal(ass2d$title, "title")
  expect_equal(ass2d$x_lab, "x_lab")
  expect_equal(ass2d$y_lab, "y_lab")
  expect_equal(ass2d$xlim, c(0, .5))
  expect_equal(ass2d$ylim, c(0, .6))
  expect_equal(ass2d$smooth_method, "loess")
  expect_equal(ass2d$smooth_benchmark, TRUE)
  expect_equal(ass2d$smooth_se, FALSE)
  expect_equal(ass2d$scale_x, "log10")
  expect_equal(ass2d$scale_y, "log10")
  expect_equal(ass2d$dpi, 250)
  expect_equal(ass2d$units, "cm")
  expect_equal(ass2d$text_size, 2.2)
  expect_equal(ass2d$alpha, 1)
  expect_equal(ass2d$width, 5)
  expect_equal(ass2d$height, 5)
  expect_equal(ass2d$colors, c("red", "blue"))
})

test_that("Ass2d$assess() works", {
  
  set.seed(132)

  n <- 25
  n_fold <- 3
  lambda <- 1
  n_genes <- 5
  n_na_in_pheno <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  model_1 <- Model$new(
    name = "cox",
    directory = file.path(dir, "cox"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    response_type = "survival_censored"
  )
  model_2 <- Model$new(
    name = "logistic",
    directory = file.path(dir, "logistic"),
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "binomial", alpha = 1, 
      nFold = n_fold, lambda = lambda, zeroSum = FALSE),
    response_type = "binary"
  )
  training_camp(list(model_1, model_2), data, quiet = TRUE)

  # ROCR metrics
  ass2d <- Ass2d$new(
    file = file.path(dir, "rpp.jpeg"),
    x_metric = "rpp",
    y_metric = "prec",
    pivot_time_cutoff = 2.,
    benchmark = "ipi",
    smooth_se = TRUE,
    show_plots = FALSE,
    text = list(ggplot2::aes(x = .5, y = .5, label = "this text"), 
     color = "red", angle = 45),
    dpi = 250,
    theme = ggplot2::theme(
      text = ggplot2::element_text(family = "Courier")
    )
  )
  ass2d$assess(data, model_1, quiet = TRUE)
  expect_s3_class(ass2d$data, "tbl_df")

  # Logrank
  ass2d$benchmark <- "ipi"
  ass2d$file <- file.path(dir, "logrank.jpeg")
  ass2d$y_metric <- "logrank"
  ass2d$scale_y <- "log10"
  ass2d$assess(data, model_2, quiet = TRUE)
  expect_equal(names(ass2d$data), c("rpp", "logrank", "cutoff", "split", "model"))

  # Precision CI
  ass2d$y_metric <- "precision_ci"
  ass2d$ci_level <- .95
  ass2d$file <- file.path(dir, "precision_ci.jpeg")
  ass2d$scale_y <- "identity"
  ass2d$title <- "Lower precision CI boundary (upper for ipi)"
  ass2d$show_plots <- FALSE
  ass2d$assess(data, model_1, quiet = TRUE)
  expect_equal(names(ass2d$data), c("rpp", "precision_ci", "cutoff", "split", "model"))
})

test_that("Ass2d$assess_center() works", {

  set.seed(4352)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  dir <- withr::local_tempdir()
  data_dir <- file.path(dir, "data")
  dir.create(data_dir, recursive = TRUE)
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = data_dir
  )
  data <- Data$new(
    name = "Mock et al. (2023)", 
    directory = data_dir, 
    train_prop = .7,
    benchmark_col = "ipi"
  )

  model_dir <- file.path(dir, "models")
  dir.create(model_dir)
  model_1 <- Model$new(
    name = "cox",
    directory = file.path(model_dir, "cox"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    response_type = "survival_censored"
  )
  model_2 <- Model$new(
    name = "binomial",
    directory = file.path(model_dir, "logistic"),
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var"
  )
  model_list <- list(model_1, model_2)

  training_camp(
    data = data,
    model_list = model_list,
    quiet = TRUE
  )

  data$cohort <- c("train", "test")
  res_dir <- file.path(dir, "results")
  ass2d <- Ass2d$new(
    file = file.path(model_dir, "perf_plot.jpeg"),
    x_metric = "rpp",
    y_metric = "prec",
    show_plots = FALSE,
    smooth_method = "loess",
    pivot_time_cutoff = 2.,
    benchmark = "ipi",
    fellow_csv = TRUE,
    scores_plot = TRUE,
    text = list(ggplot2::aes(x = .5, y = .5, label = "hello"), 
      color = "red", angle = 90),
    theme = ggplot2::theme_minimal() + 
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "red")),
    dpi = 250
  )
  ass2d$assess_center(data, model_list, comparison_plot = TRUE, risk_scores = TRUE, 
    quiet = TRUE)
  expect_true(file.exists(ass2d$file))
  expect_true(file.exists(file.path(res_dir, "logistic/2/scores.jpeg")))
  expect_true(file.exists(file.path(res_dir, "cox/2/rpp_vs_prec.csv")))
  expect_true(file.exists(file.path(model_dir, "perf_plot.jpeg")))

  data$cohort <- "test"
  model_1$split_index <- 1
  ass2d$y_metric <- "logrank"
  ass2d$scale_y <- "log10"
  ass2d$file <- file.path(model_dir, "logrank.jpeg")
  ass2d$benchmark <- NULL
  ass2d$text <- NULL
  ass2d$scores_plot <- FALSE
  ass2d$assess_center(data, model_list, comparison_plot = TRUE, quiet = TRUE)
  expect_true(file.exists(file.path(res_dir, "logrank.jpeg")))

  data$cohort <- "train"
  model_2$time_cutoffs <- 1.5
  ass2d$y_metric <- "precision_ci"
  ass2d$ci_level <- .95
  ass2d$file <- file.path(model_dir, "precision_ci.jpeg")
  ass2d$benchmark <- "ipi"
  ass2d$assess_center(data, model_list, comparison_plot = FALSE, quiet = TRUE)
  expect_true(file.exists(file.path(model_dir, "logistic/1-5/rpp_vs_precision_ci.jpeg")))
})
