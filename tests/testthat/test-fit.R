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
  save_model_as <- c("mock", "this_model")
  additional_model_par <- list("alpha" = 0.5, nFold = 3)
  results_dir <- withr::local_tempdir(pattern = "results")
  
  # Case 1:
  model <- "cox_lasso_zerosum"

  # non-existing directory
  expect_error(fit(
    x, y, model, save_model_as, 
    additional_model_par, 
    results_dir = results_dir , 
    create_dir = FALSE
    ),
    regexp = "not exist")
  
  # valid
  expect_message(fit_obj <- fit(
    x, y, model, save_model_as, 
    additional_model_par = additional_model_par, 
    results_dir = results_dir,
    create_dir = TRUE,
    plots = TRUE
    ),
    regexp = "Creating directory")
  expect_s3_class(fit_obj, "zeroSum")
  
  # Case 2:
  model <- "invalid_model" 

  expect_error(fit(
    x, y, "invalid_model", save_model_as, 
    additional_model_par = additional_model_par, 
    results_dir = results_dir
    ),
    regexp = "not supported")
  
  # Case 3: model = "lasso_zerosum"
  y <- y[, 2, drop = FALSE]
  model <- "lasso_zerosum"
  expect_message(fit_obj <- fit(
      x, y, model, save_model_as, 
      additional_model_par = additional_model_par, 
      results_dir = results_dir,
      create_dir = TRUE,
      plots = FALSE
      ),
      regexp = "Creating directory")
  expect_s3_class(fit_obj, "zeroSum")
})