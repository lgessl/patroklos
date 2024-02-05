test_that("calculate_2d_metric() works", {
    
  set.seed(134)

  n_samples <- 10
  split_index <- 1:2

  l <- apb(n_samples, split_index, fluctuating_availability = FALSE)
  actual <- l[[1]]
  predicted <- l[[2]]
  benchmark <- l[[3]]
  
  model_spec <- ModelSpec(
    name = "mock model",
    directory = "mock",
    fitter = zeroSum::zeroSum,
    split_index = split_index,
    time_cutoffs = 2.
  )
  ass_2d_spec <- Ass2dSpec(
    fname = "test.pdf",
    x_metric = "rpp",
    y_metric = "prec",
    benchmark = "ipi"
  )
  pheno_tbl <- generate_mock_data(
    n_samples = n_samples,
    n_genes = 2,
    n_na_in_pheno = 0
  )[["pheno_tbl"]]
  data_spec <- DataSpec(
    name = "mock data",
    directory = "mock",
    train_prop = .66
  )

  ass_2d_spec <- calculate_2d_metric(
    actual = actual,
    predicted = predicted,
    benchmark = benchmark,
    ass_2d_spec = ass_2d_spec,
    model_spec = model_spec
  )
  perf_tbl <- ass_2d_spec$data
  expect_equal(names(perf_tbl), c("rpp", "prec", "cutoff", "split", "model"))
  expect_s3_class(perf_tbl, "tbl_df")
  expect_true(all(perf_tbl[["model"]] %in% c(model_spec$name, ass_2d_spec$benchmark)))
  expect_true(perf_tbl[, 1:4] |> as.matrix() |> is.numeric() |> all())

  ass_2d_spec$y_metric <- "logrank"
  ass_2d_spec <- calculate_2d_metric(
    actual = actual,
    predicted = predicted,
    ass_2d_spec = ass_2d_spec,
    model_spec = model_spec,
    benchmark = benchmark,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec
  )

  ass_2d_spec$y_metric <- "precision_ci"
  ass_2d_spec$ci_level <- .95
  ass_2d_spec <- calculate_2d_metric(
    actual = actual,
    predicted = predicted,
    ass_2d_spec = ass_2d_spec,
    model_spec = model_spec,
    benchmark = benchmark,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec
  )
})


test_that("precision_ci() works", {

  set.seed(143)

  n_samples <- 10
  n_extra <- 2
  gamma <- .85

  estimate <- rnorm(n_samples)
  actual <- sample(c(0, 1), n_samples + n_extra, replace = TRUE)
  names(estimate) <- paste0("s", sample(n_samples))
  names(actual) <- paste0("s", sample(length(actual)))
  actual[1] <- NA

  tbl_low <- precision_ci(
    estimate = estimate,
    actual = actual,
    confidence_level = gamma,
    y_metric = "ci_boundary",
    x_metric = "prevalence",
    lower_boundary = TRUE
  )
  tbl_high <- precision_ci(
    estimate = estimate,
    actual = actual,
    confidence_level = gamma,
    y_metric = "ci_boundary",
    x_metric = "prevalence",
    lower_boundary = FALSE
  )
  expect_true(all(tbl_low[, 1] <= tbl_high[, 1]))
})


test_that("metric_with_rocr() works", {

  set.seed(354)

  n_samples <- 10
  n_extra <- 1

  estimate <- rnorm(n_samples)
  actual <- sample(c(0, 1), n_samples + n_extra, replace = TRUE)
  names(estimate) <- paste0("s", sample(n_samples))
  names(actual) <- paste0("s", sample(length(actual)))
  actual[1] <- NA

  tbl <- metric_with_rocr(
    estimate = estimate,
    actual = actual,
    x_metric = "rpp",
    y_metric = "prec"
  )
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), 3)
})


test_that("logrank_metric() works", {

  set.seed(354)

  n_samples <- 10
  n_extra <- 5
  n_na_in_estimate <- 1

  pheno_tbl <- generate_mock_data(
    n_samples = n_samples + n_extra,
    n_genes = 2,
    n_na_in_pheno = 0
  )[["pheno_tbl"]]
  data_spec <- DataSpec(
    name = "mock",
    directory = "mock",
    train_prop = .66
  )

  estimate <- rnorm(n_samples)
  estimate[sample(length(estimate), n_na_in_estimate)] <- NA
  names(estimate) <- sample(pheno_tbl[["patient_id"]], n_samples)

  tbl <- logrank_metric(
    estimate = estimate,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    y_metric = "my_y",
    x_metric = "my_x"
  )
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), 3)
  expect_equal(nrow(tbl), length(estimate) - 1 - n_na_in_estimate)
})
