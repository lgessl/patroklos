test_that("hypertune() and its return value's predict method work", {

  set.seed(123)

  n_samples <- 10
  n_genes <- 2

  x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")

  tune_fitter <- hypertune(ptk_ranger)
  ptk_hypertune <- tune_fitter(x3y[[1]], x3y[2:4], val_error_fun = error_rate, 
    num.trees = c(30, 35), rel_mtry = FALSE, mtry = c(2, 3), classification = TRUE)
  expect_s3_class(ptk_hypertune, "ptk_hypertune")
  expect_true(is.character(ptk_hypertune$lambda) && length(ptk_hypertune$lambda) == 4)
  expect_true(is.numeric(ptk_hypertune$val_error) && length(ptk_hypertune$val_error) == 4)
  expect_true(is.numeric(ptk_hypertune$lambda_min_index) && length(ptk_hypertune$lambda_min_index) == 1)
  expect_true(is.character(ptk_hypertune$lambda_min) && length(ptk_hypertune$lambda_min) == 1)
  expect_true(is.numeric(ptk_hypertune$min_error) && length(ptk_hypertune$min_error) == 1)
  expect_true(is.list(ptk_hypertune$fit_obj_list) && length(ptk_hypertune$fit_obj_list) == 4)
  expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) all(x %in% c(0, 1)))))
  expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) inherits(x, "ranger"))))
  expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) nrow(x) == sum(!is.na(x3y[[2]])))))
  expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) all(dim(x) == c(n_samples, 1)))))
  expect_equal(ptk_hypertune$val_predict_list[[ptk_hypertune$lambda_min_index]], 
    ptk_hypertune$val_predict)

  pred <- predict(ptk_hypertune, newx = x3y[[1]])
  expect_equal(length(pred), n_samples)

  # tune_fitter <- hypertune(ptk_zerosum)
  # ptk_hypertune <- tune_fitter(x3y[[1]], x3y[2:4], val_error_fun = neg_binomial_log_likelihood, 
  #   lambda = c(0.5, 1), nFold = 2, family = "binomial")
  # expect_s3_class(ptk_hypertune, "ptk_hypertune")
  # expect_true(is.character(ptk_hypertune$lambda) && length(ptk_hypertune$lambda) == 2)
  # expect_true(is.numeric(ptk_hypertune$val_error) && length(ptk_hypertune$val_error) == 2)
  # expect_true(is.numeric(ptk_hypertune$lambda_min_index) && length(ptk_hypertune$lambda_min_index) == 1)
  # expect_true(is.character(ptk_hypertune$lambda_min) && length(ptk_hypertune$lambda_min) == 1)
  # expect_true(is.numeric(ptk_hypertune$min_error) && length(ptk_hypertune$min_error) == 1)
  # expect_true(is.character(ptk_hypertune$error_name) && length(ptk_hypertune$error_name) == 1)
  # expect_true(is.list(ptk_hypertune$fit_obj_list) && length(ptk_hypertune$fit_obj_list) == 2)
  # expect_true(all(sapply(ptk_hypertune$val_predict_list, function(x) nrow(x) == sum(!is.na(x3y[[2]])))))
  # expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) all(dim(x) == c(n_samples, 1)))))
  # expect_true(all(sapply(ptk_hypertune$fit_obj_list, function(x) inherits(x, "ptk_zerosum"))))

  # pred <- predict(ptk_hypertune, newx = x3y[[1]])
  # expect_equal(length(pred), n_samples)

  # Handling of NA returned by fitter
  tune_fitter <- hypertune(ptk_ranger)
  ptk_hypertune <- tune_fitter(x3y[[1]], x3y[2:4], val_error_fun = error_rate, num.trees = c(30, 35), 
    rel_mtry = FALSE, mtry = c(2, 20), classification = TRUE, skip_on_invalid_input = TRUE)
  expect_equal(length(ptk_hypertune$fit_obj_list), 2)

  # select = TRUE
  tune_fitter <- hypertune(ptk_ranger, select = TRUE)
  ptk_hypertune <-  tune_fitter(x3y[[1]], x3y[2:4], val_error_fun = error_rate, 
    num.trees = c(30, 35), rel_mtry = FALSE, mtry = c(2, 3), classification = TRUE)
  expect_true(all(is.na(ptk_hypertune$fit_obj_list[-ptk_hypertune$lambda_min_index])))
  expect_true(is.na(ptk_hypertune$val_predict_list))
  expect_s3_class(ptk_hypertune$fit_obj_list[[ptk_hypertune$lambda_min_index]], "ranger")
})

test_that("unitune() works", {

  set.seed(213)

  n_samples <- 10
  n_genes <- 2

  x3y <- generate_mock_data(n_samples = n_samples, n_genes = n_genes, 
    return_type = "fitter")
  combos <- matrix(sample(0:1, n_samples * 3, replace = TRUE), ncol = 3)
  rownames(combos) <- rownames(x3y[[1]])
  colnames(combos) <- c("a&b", "a&c", "a&b&c")
  x3y[[1]] <- cbind(x3y[[1]], combos)
  attr(x3y[[1]], "li_var_suffix") <- "++"

  fit_obj <- unitune(ptk_zerosum)(x3y[[1]], x3y[[3]], time_cutoffs = c(1.75, 2),
    val_error_fun = neg_roc_auc, pivot_time_cutoff = 2, 
    combine_n_max_categorical_features = 2:3, lambdaSteps = 2, nFold = 2, 
    family = "binomial", zeroSum = FALSE)
  expect_s3_class(fit_obj, "ptk_zerosum")
  expect_true(fit_obj$combine_n_max_categorical_features %in% c(2, 3))
  expect_equal(length(fit_obj$val_predict_list), 2)
  expect_equal(length(unique(fit_obj$foldid)), 2)
  expect_true("a&b" %in% fit_obj$variables.names)

  fit_obj <- unitune(ptk_zerosum)(x3y[[1]], x3y[[3]], time_cutoffs = c(2, 2.25),
    val_error_fun = neg_binomial_log_likelihood, pivot_time_cutoff = 2, 
    combine_n_max_categorical_features = 1:2, lambdaSteps = 2, nFold = 2, 
    family = "binomial", zeroSum = FALSE)
  expect_false("a&b&c" %in% fit_obj$variables.names)
})
