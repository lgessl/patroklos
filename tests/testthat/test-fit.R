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
  optional_fitter_args <- list("alpha" = 0.5, nFold = 3)
  results_dir <- withr::local_tempdir(pattern = "results")
  save_dir_base <- file.path(results_dir, "mock")
  
  fitter <- zeroSum::zeroSum

  # Case 1:
  save_dir <- file.path(save_dir_base, "cox_lasso_zerosum", "std")
  optional_fitter_args[["family"]] <- "cox"

  # non-existing directory
  expect_error(fit(
      x, y, fitter, save_dir, 
      optional_fitter_args,  
      create_save_dir = FALSE
    ),
    regexp = "not exist"
  )
  
  # valid
  expect_message(fit_obj <- fit(
      x, y, fitter, save_dir, 
      optional_fitter_args = optional_fitter_args, 
      create_save_dir = TRUE,
      save_plots = TRUE
    ),
    regexp = "Creating directory"
  )
  expect_s3_class(fit_obj, "zeroSum")
  
  # Case 2:
  fake_fitter <- function(x) {}

  expect_error(fit(
      x, y, fake_fitter, save_dir, 
      optional_fitter_args = optional_fitter_args
    )
  )
  
  # Case 3: model = "lasso_zerosum"
  save_dir <- file.path(save_dir_base, "lasso_zerosum", "std")
  optional_fitter_args[["family"]] <- "binomial"
  y <- y[, 2, drop = FALSE]

  expect_message(fit_obj <- fit(
      x, y, fitter, save_dir, 
      optional_fitter_args = optional_fitter_args, 
      create_save_dir = TRUE,
      save_plots = FALSE,
      plot_fname = "other.pdf",
      fit_fname = "other.rds"
    ),
    regexp = "Creating directory"
  )
  expect_s3_class(fit_obj, "zeroSum")
})