test_that("hypertune() and its return value's predict method work", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2

  xyy <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")

  tune_fitter <- hypertune(ptk_ranger, error = "error_rate")
  ptk_hypertune <- tune_fitter(xyy[[1]], xyy[[2]], xyy[[3]], num.trees = c(30, 35), 
    rel_mtry = FALSE, mtry = c(2, 3), classification = TRUE)
  expect_s3_class(ptk_hypertune, "ptk_hypertune")
  expect_true(is.character(ptk_hypertune$lambda) && length(ptk_hypertune$lambda) == 4)
  expect_true(is.numeric(ptk_hypertune$val_error) && length(ptk_hypertune$val_error) == 4)
  expect_true(is.numeric(ptk_hypertune$lambda_min_index) && length(ptk_hypertune$lambda_min_index) == 1)
  expect_true(is.character(ptk_hypertune$min_lambda) && length(ptk_hypertune$min_lambda) == 1)
  expect_true(is.numeric(ptk_hypertune$min_error) && length(ptk_hypertune$min_error) == 1)
  expect_true(is.character(ptk_hypertune$error_name) && length(ptk_hypertune$error_name) == 1)
  expect_true(is.list(ptk_hypertune$fit_obj_list) && length(ptk_hypertune$fit_obj_list) == 4)
  expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) all(x %in% c(0, 1)))))
  expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) inherits(x, "ranger"))))
  expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) nrow(x) == sum(!is.na(xyy[[2]])))))
  expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) all(dim(x) == c(n_samples, 1)))))

  pred <- predict(ptk_hypertune, newx = xyy[[1]])
  expect_equal(length(pred), n_samples)

  tune_fitter <- hypertune(ptk_zerosum, error = "neg_roc_auc")
  ptk_hypertune <- tune_fitter(xyy[[1]], xyy[[2]], xyy[[3]], lambda = c(0.5, 1), 
    nFold = 2, family = "binomial")
  expect_s3_class(ptk_hypertune, "ptk_hypertune")
  expect_true(is.character(ptk_hypertune$lambda) && length(ptk_hypertune$lambda) == 2)
  expect_true(is.numeric(ptk_hypertune$val_error) && length(ptk_hypertune$val_error) == 2)
  expect_true(is.numeric(ptk_hypertune$lambda_min_index) && length(ptk_hypertune$lambda_min_index) == 1)
  expect_true(is.character(ptk_hypertune$min_lambda) && length(ptk_hypertune$min_lambda) == 1)
  expect_true(is.numeric(ptk_hypertune$min_error) && length(ptk_hypertune$min_error) == 1)
  expect_true(is.character(ptk_hypertune$error_name) && length(ptk_hypertune$error_name) == 1)
  expect_true(is.list(ptk_hypertune$fit_obj_list) && length(ptk_hypertune$fit_obj_list) == 2)
  expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) nrow(x) == sum(!is.na(xyy[[2]])))))
  expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) all(dim(x) == c(n_samples, 1)))))
  expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) inherits(x, "ptk_zerosum"))))

  pred <- predict(ptk_hypertune, newx = xyy[[1]])
  expect_equal(length(pred), n_samples)

  # Handling of NA returned by fitter
  tune_fitter <- hypertune(ptk_ranger, error = "neg_roc_auc")
  ptk_hypertune <- tune_fitter(xyy[[1]], xyy[[2]], xyy[[3]], num.trees = c(30, 35), 
    rel_mtry = FALSE, mtry = c(2, 20), classification = TRUE, 
    skip_on_invalid_input = TRUE)
  expect_equal(length(ptk_hypertune$fit_obj_list), 2)
})
