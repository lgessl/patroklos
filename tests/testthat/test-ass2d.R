test_that("Ass2d$new() works", {

  ass2d <- Ass2d$new(
    x_metric = "rpp",
    y_metric = "prec",
    ci_level = .66,
    x_lab = "x_lab",
    y_lab = "y_lab",
    xlim = c(0, .5),
    ylim = c(0, .6),
    scale_x = "log10",
    scale_y = "log10",
    dpi = 250,
    units = "cm",
    alpha = 1,
    width = 5,
    height = 5,
    colors = c("red", "blue")
  )
  expect_s3_class(ass2d, "Ass2d")
  expect_equal(ass2d$x_metric, "rpp")
  expect_equal(ass2d$y_metric, "prec")
  expect_equal(ass2d$ci_level, .66)
  expect_equal(ass2d$x_lab, "x_lab")
  expect_equal(ass2d$y_lab, "y_lab")
  expect_equal(ass2d$xlim, c(0, .5))
  expect_equal(ass2d$ylim, c(0, .6))
  expect_equal(ass2d$scale_x, "log10")
  expect_equal(ass2d$scale_y, "log10")
  expect_equal(ass2d$dpi, 250)
  expect_equal(ass2d$units, "cm")
  expect_equal(ass2d$alpha, 1)
  expect_equal(ass2d$width, 5)
  expect_equal(ass2d$height, 5)
  expect_equal(ass2d$colors, c("red", "blue"))
})

test_that("Ass2d$assess() works", {
  
  set.seed(132)

  n <- 40
  n_fold <- 3
  n_genes <- 5
  n_na_in_pheno <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  model <- Model$new(
    name = "logistic",
    directory = file.path(dir, "logistic"),
    fitter = ptk_zerosum,
    time_cutoffs = 2.,
    val_error_fun = neg_prec_with_prev_greater(0.15),
    hyperparams = list(family = "binomial", alpha = 1, 
      nFold = n_fold, lambda = 0, zeroSum = FALSE)
  )
  training_camp(list(model), data, quiet = TRUE, skip_on_error = FALSE)

  # ROCR metrics
  ass2d <- Ass2d$new(
    x_metric = "rpp",
    y_metric = "prec",
    dpi = 250,
    theme = ggplot2::theme(
      text = ggplot2::element_text(family = "Courier")
    )
  )
  expect_no_error(
    plt <- ass2d$assess(data, model, quiet = TRUE)
  )
  # print(plt)

  # Logrank
  ass2d$y_metric <- "logrank"
  ass2d$y_lab <- "logrank p"
  ass2d$scale_y <- "log10"
  expect_no_error(
    plt <- ass2d$assess(data, model, quiet = TRUE)
  )
  # print(plt)

  # Precision CI
  ass2d$y_metric <- "precision_ci"
  ass2d$ci_level <- .95
  ass2d$scale_y <- "identity"
  expect_no_error(
    plt <- ass2d$assess(data, model, quiet = TRUE)
  )
  # print(plt)
})

test_that("Ass2d$assess_center() works", {

  set.seed(4352)

  n_samples <- 50
  n_genes <- 2
  n_na_in_pheno <- 5

  dir <- withr::local_tempdir()
  data_dir <- file.path(dir, "data")
  dir.create(data_dir, recursive = TRUE)
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = data_dir
  )

  model_dir <- file.path(dir, "models")
  dir.create(model_dir)
  model_1 <- Model$new(
    name = "cox",
    directory = file.path(model_dir, "cox"),
    fitter = ptk_zerosum,
    time_cutoffs = Inf,
    val_error_fun = neg_roc_auc,
    hyperparams = list(family = "cox", alpha = 1, nFold = 2, 
      lambda = 0, zeroSum = FALSE)
  )
  model_2 <- Model$new(
    name = "binomial",
    directory = file.path(model_dir, "logistic"),
    fitter = ptk_zerosum,
    hyperparams = list(family = "binomial", alpha = 1, nFold = 2, 
      lambda = 0, zeroSum = FALSE),
    time_cutoffs = c(1.5, 2),
    val_error_fun = neg_binomial_log_likelihood,
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var"
  )
  model_3 <- Model$new(
    name = "rf",
    directory = file.path(model_dir, "rf"),
    fitter = multitune(ptk_ranger),
    hyperparams = list(rel_mtry = FALSE, mtry = 2, num.trees = 20, 
      min.node.size = 5, classification = TRUE),
    time_cutoffs = 2,
    val_error_fun = error_rate,
    include_from_discrete_pheno = c("discrete_var", "abc_gcb", "ipi_age"),
    combine_n_max_categorical_features = 3,
    combined_feature_min_positive_ratio = 0.1
  )
  model_list <- list(model_1, model_2, model_3)

  training_camp(
    data = data,
    model_list = model_list,
    quiet = TRUE,
    skip_on_error = FALSE
  )

  res_dir <- file.path(dir, "results")
  ass2d <- Ass2d$new(
    x_metric = "rpp",
    y_metric = "prec",
    theme = ggplot2::theme_minimal() + 
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")),
    dpi = 250
  )
  data$cohort <- "test"
  expect_no_error(
    plt <- ass2d$assess_center(data, model_list, quiet = TRUE)
  )
  # print(plt)

  ass2d$y_metric <- "logrank"
  ass2d$scale_y <- "log10"
  expect_no_error(
    plt <- ass2d$assess_center(data, model_list, quiet = TRUE)
  )
  # print(plt)

  data$cohort <- "val_predict"
  ass2d$x_metric <- "rank"
  ass2d$y_metric <- "risk score"
  ass2d$ci_level <- .95
  ass2d$scale_y <- "identity"
  expect_no_error(
    plt <- ass2d$assess_center(data, model_list, quiet = TRUE)
  )
  # print(plt)
})
