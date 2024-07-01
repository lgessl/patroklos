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
  data <- Data$new(
    name = "mock",
    directory = "mock",
    pivot_time_cutoff = 2,
    time_to_event_col = "pfs_years",
    cohort = "train"
  )

  estimate <- rnorm(n_samples)
  estimate[sample(length(estimate), n_na_in_estimate)] <- NA
  names(estimate) <- sample(pheno_tbl[["patient_id"]], n_samples)

  tbl <- logrank_metric(
    estimate = estimate,
    pheno_tbl = pheno_tbl,
    data = data,
    y_metric = "my_y",
    x_metric = "my_x"
  )
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), 3)
  expect_equal(nrow(tbl), length(estimate) - 1 - n_na_in_estimate)
})
