test_that("precision_ci() works", {

  set.seed(143)

  n_samples <- 10
  n_extra <- 2
  gamma <- .85

  predicted <- rnorm(n_samples)
  actual <- sample(c(0, 1), n_samples, replace = TRUE)

  tbl_low <- precision_ci(
    predicted = predicted,
    actual = actual,
    confidence_level = gamma
  )
  tbl_high <- precision_ci(
    predicted = predicted,
    actual = actual,
    confidence_level = gamma
  )
  expect_true(all(tbl_low[, 1] <= tbl_high[, 1]))
})


test_that("metric_with_rocr() works", {

  set.seed(354)

  n_samples <- 10

  predicted <- rnorm(n_samples)
  actual <- sample(c(0, 1), n_samples, replace = TRUE)
  names(predicted) <- paste0("s", sample(n_samples))
  names(actual) <- names(predicted)

  tbl <- metric_with_rocr(
    predicted = predicted,
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

  predicted <- rnorm(n_samples)
  names(predicted) <- sample(pheno_tbl[["patient_id"]], n_samples)

  tbl <- logrank_metric(
    predicted = predicted,
    pheno_tbl = pheno_tbl,
    data = data
  )
  expect_s3_class(tbl, "tbl_df")
  expect_equal(ncol(tbl), 3)
  expect_equal(nrow(tbl), length(predicted)-1)
})