test_that("integrated_cv() and its return value's predict method work", {

    set.seed(123)

    n_samples <- 10
    n_genes <- 2

    x_y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
      return_type = "fitter")
    x <- x_y$x
    y <- x_y$y

    cv_fitter <- integrated_cv(ptk_ranger, metric = "accuracy")
    ptk_cv <- cv_fitter(x, y, num.trees = c(30, 35), mtry = c(2, 3), 
      classification = TRUE)
    expect_s3_class(ptk_cv, "ptk_cv")
    expect_true(is.character(ptk_cv$lambda) && length(ptk_cv$lambda) == 4)
    expect_true(is.numeric(ptk_cv$val_metric) && length(ptk_cv$val_metric) == 4)
    expect_true(is.numeric(ptk_cv$best_lambda_index) && length(ptk_cv$best_lambda_index) == 1)
    expect_true(is.character(ptk_cv$best_lambda) && length(ptk_cv$best_lambda) == 1)
    expect_true(is.numeric(ptk_cv$best_metric) && length(ptk_cv$best_metric) == 1)
    expect_true(is.character(ptk_cv$metric_name) && length(ptk_cv$metric_name) == 1)
    expect_true(is.list(ptk_cv$fit_obj_list) && length(ptk_cv$fit_obj_list) == 4)
    expect_true(all(sapply(ptk_cv$val_predict_list, function(x) all(x %in% c(0, 1)))))
    expect_true(all(sapply(ptk_cv$fit_obj_list, function(x) inherits(x, "ranger"))))

    cv_fitter <- integrated_cv(ptk_zerosum, metric = "roc_auc")
    ptk_cv <- cv_fitter(x, y, lambda = c(0.5, 1), nFold = 2, family = "binomial")
    expect_s3_class(ptk_cv, "ptk_cv")
    expect_true(is.character(ptk_cv$lambda) && length(ptk_cv$lambda) == 2)
    expect_true(is.numeric(ptk_cv$val_metric) && length(ptk_cv$val_metric) == 2)
    expect_true(is.numeric(ptk_cv$best_lambda_index) && length(ptk_cv$best_lambda_index) == 1)
    expect_true(is.character(ptk_cv$best_lambda) && length(ptk_cv$best_lambda) == 1)
    expect_true(is.numeric(ptk_cv$best_metric) && length(ptk_cv$best_metric) == 1)
    expect_true(is.character(ptk_cv$metric_name) && length(ptk_cv$metric_name) == 1)
    expect_true(is.list(ptk_cv$fit_obj_list) && length(ptk_cv$fit_obj_list) == 2)
    expect_true(all(sapply(ptk_cv$val_predict_list, function(x) length(x) == n_samples)))
    expect_true(all(sapply(ptk_cv$fit_obj_list, function(x) inherits(x, "ptk_zerosum"))))

    pred <- predict(ptk_cv, newx = x)
    expect_equal(length(pred), n_samples)
})