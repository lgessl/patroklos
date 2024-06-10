test_that("hypertune() and its return value's predict method work", {

    set.seed(123)

    n_samples <- 10
    n_genes <- 2

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
      return_type = "fitter")
    x <- x_y$x
    y <- x_y$y

    tune_fitter <- hypertune(ptk_ranger, metric = "accuracy")
    ptk_hypertune <- tune_fitter(x, y, num.trees = c(30, 35), mtry = c(2, 3), 
      classification = TRUE)
    expect_s3_class(ptk_hypertune, "ptk_hypertune")
    expect_true(is.character(ptk_hypertune$lambda) && length(ptk_hypertune$lambda) == 4)
    expect_true(is.numeric(ptk_hypertune$val_metric) && length(ptk_hypertune$val_metric) == 4)
    expect_true(is.numeric(ptk_hypertune$best_lambda_index) && length(ptk_hypertune$best_lambda_index) == 1)
    expect_true(is.character(ptk_hypertune$best_lambda) && length(ptk_hypertune$best_lambda) == 1)
    expect_true(is.numeric(ptk_hypertune$best_metric) && length(ptk_hypertune$best_metric) == 1)
    expect_true(is.character(ptk_hypertune$metric_name) && length(ptk_hypertune$metric_name) == 1)
    expect_true(is.list(ptk_hypertune$fit_obj_list) && length(ptk_hypertune$fit_obj_list) == 4)
    expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) all(x %in% c(0, 1)))))
    expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) inherits(x, "ranger"))))

    pred <- predict(ptk_hypertune, newx = x)
    expect_equal(length(pred), n_samples)

    tune_fitter <- hypertune(ptk_zerosum, metric = "roc_auc")
    ptk_hypertune <- tune_fitter(x, y, lambda = c(0.5, 1), nFold = 2, family = "binomial")
    expect_s3_class(ptk_hypertune, "ptk_hypertune")
    expect_true(is.character(ptk_hypertune$lambda) && length(ptk_hypertune$lambda) == 2)
    expect_true(is.numeric(ptk_hypertune$val_metric) && length(ptk_hypertune$val_metric) == 2)
    expect_true(is.numeric(ptk_hypertune$best_lambda_index) && length(ptk_hypertune$best_lambda_index) == 1)
    expect_true(is.character(ptk_hypertune$best_lambda) && length(ptk_hypertune$best_lambda) == 1)
    expect_true(is.numeric(ptk_hypertune$best_metric) && length(ptk_hypertune$best_metric) == 1)
    expect_true(is.character(ptk_hypertune$metric_name) && length(ptk_hypertune$metric_name) == 1)
    expect_true(is.list(ptk_hypertune$fit_obj_list) && length(ptk_hypertune$fit_obj_list) == 2)
    expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) length(x) == n_samples)))
    expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) inherits(x, "ptk_zerosum"))))

    pred <- predict(ptk_hypertune, newx = x)
    expect_equal(length(pred), n_samples)

  # Handling of NA returned by fitter
  tune_fitter <- hypertune(ptk_ranger, metric = "accuracy")
  ptk_hypertune <- tune_fitter(x, y, num.trees = c(30, 35), mtry = c(2, 20), 
    classification = TRUE, skip_on_invalid_input = TRUE)
  expect_equal(length(ptk_hypertune$fit_obj_list), 2)
})
