test_that("assess_2d_center() works", {

  set.seed(4352)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  base_dir <- withr::local_tempdir("data")
  data_dir <- file.path(base_dir, "data/mock")
  dir.create(data_dir, recursive = TRUE)
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = data_dir
  )
  data <- Data(
    name = "Mock et al. (2023)", 
    directory = data_dir, 
    train_prop = .7,
    benchmark_col = "ipi"
  )

  model_dir <- file.path(base_dir, "models")
  dir.create(model_dir)
  model_1 <- Model$new(
    name = "cox-zerosum",
    directory = file.path(model_dir, "cox"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    response_type = "survival_censored",
    fit_file = "model1.rds"
  )
  model_2 <- Model$new(
    name = "binomial-zerosum",
    directory = file.path(model_dir, "logistic"),
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    split_index = 1,
    time_cutoffs = c(1.5, 2),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    fit_file = "model2.rds"
  )
  model_list <- list(model_1, model_2)

  training_camp(
    data = data,
    model_list = model_list,
    quiet = TRUE
  )

  res_dir <- file.path(base_dir, "results")
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

  assess_2d_center(
    model_list = model_list,
    data = data,
    ass2d = ass2d,
    cohort = c("train", "test"),
    quiet = TRUE
  )
  expect_true(file.exists(ass2d$file))
  expect_true(file.exists(file.path(res_dir, "logistic/2/scores.jpeg")))
  expect_true(file.exists(file.path(res_dir, "cox/2/rpp_vs_prec.csv")))

  model_1$split_index <- 1
  ass2d$y_metric <- "logrank"
  ass2d$scale_y <- "log10"
  ass2d$file <- file.path(model_dir, "logrank.jpeg")
  ass2d$benchmark <- NULL
  ass2d$text <- NULL
  ass2d$scores_plot <- FALSE
  assess_2d_center(
    model_list = model_list,
    data = data,
    ass2d = ass2d,
    comparison_plot = FALSE,
    cohorts = "test",
    quiet = TRUE
  )

  model_2$time_cutoffs <- 1.5
  ass2d$y_metric <- "precision_ci"
  ass2d$ci_level <- .95
  ass2d$file <- file.path(model_dir, "precision_ci.jpeg")
  ass2d$benchmark <- "ipi"
  assess_2d_center(
    model_list = list(model_2),
    data = data,
    ass2d = ass2d,
    comparison_plot = TRUE,
    cohorts = "train",
    quiet = TRUE
  )
  expect_true(file.exists(file.path(model_dir, "logistic/1-5/rpp_vs_precision_ci.jpeg")))
})
