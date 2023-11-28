library(testthat)

test_that("fit function works correctly", {

  set.seed(4325)
  n_samples = 5
  n_features = 5

  n_fold = floor(n_samples/2)
  x <- matrix(rnorm(n_samples*n_features), nrow = n_samples)
  y <- rnorm(n_samples) |> as.matrix()
  y <- cbind(y, sample(c(0, 1), n_samples, replace = TRUE))
  rownames(y) <- stringr::str_c("sample_", 1:n_samples)
  rownames(x) <- rownames(y)
  results_dir <- withr::local_tempdir(pattern = "results")
  save_dir_base <- file.path(results_dir, "mock")

  model_spec <- ModelSpec(
    fitter = zeroSum::zeroSum,
    response_type = "survival_censored",
    optional_fitter_args = list("alpha" = 0.5, nFold = 3),
    pfs_leq = 2.0,
    plot_fname = "foo.pdf",
    fit_fname = "bar.rds"
  )

  # Case 1:
  model_spec$save_dir <- file.path(save_dir_base, "cox_lasso_zerosum", "std")
  model_spec$optional_fitter_args[["family"]] <- "cox"
  model_spec$create_save_dir <- FALSE

  # non-existing directory
  expect_error(fit(x, y, model_spec), regexp = "not exist")
  
  # valid
  model_spec$create_save_dir <- TRUE
  expect_message(fit_obj <- fit(x, y, model_spec), regexp = "Creating directory")
  expect_s3_class(fit_obj, "zeroSum")
  
  # Case 2:
  model_spec$fitter <- function(x) {}

  expect_error(fit(x, y, model_spec))
  
  # Case 3: model = "lasso_zerosum"
  model_spec$fitter <- zeroSum::zeroSum
  model_spec$save_dir <- file.path(save_dir_base, "lasso_zerosum", "std")
  model_spec$optional_fitter_args[["family"]] <- "binomial"
  y <- y[, 2, drop = FALSE]

  expect_message(fit_obj <- fit(x, y, model_spec), regexp = "Creating directory")
  expect_s3_class(fit_obj, "zeroSum")
})