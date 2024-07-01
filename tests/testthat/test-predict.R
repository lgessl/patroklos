test_that("model$predict() works", {

  set.seed(4352)

  n_samples <- 50
  n_genes <- 5
  n_na_in_pheno <- 5
  n_fold <- 1
  lambda <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )

  model1 <- Model$new(
    name = "cox-zerosum",
    directory = dir,
    fitter = ptk_zerosum,
    time_cutoffs = Inf,
    val_error_fun = neg_roc_auc,
    hyperparams = list(family = "cox", alpha = 1, nFold = n_fold, 
    lambda = lambda, zeroSum = FALSE)
  )
  model1$fit(data, quiet = TRUE)
  model2 <- Model$new(
    name = "log-rf",
    directory = dir,
    fitter = long_nestor,
    time_cutoffs = 2,
    val_error_fun = neg_prec_with_prev_greater(0.15),
    hyperparams = list(
      fitter1 = ptk_zerosum,
      fitter2 = hypertune(ptk_ranger),
      hyperparams1 = list(family = "binomial", alpha = 1, nFold = n_fold, 
        lambdaSteps = 1, zeroSum = FALSE),
      hyperparams2 = list(rel_mtry = FALSE, mtry = 3, min.node.size = 5,
        classification = TRUE, num.trees = 100, skip_on_invalid_input = TRUE)
    )
  )
  model2$fit(data, quiet = TRUE)

  data$cohort <- "test"
  res <- model1$predict(data, quiet = TRUE)
  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(sapply(res, is.numeric)))
  expect_true(all(res[["actual"]] %in% c(0, 1, NA)))
  expect_true(all(sapply(res, is.vector)))
  expect_true(all(sapply(res, function(x) !is.null(names(x)))))

  data$cohort <- "val_predict"
  res <- model2$predict(data, quiet = TRUE)
  expect_true(is.list(res))
  expect_equal(length(res), 3)
  expect_true(all(sapply(res[1:2], is.numeric))) # third list element is NULL
  expect_true(all(res[["actual"]] %in% c(0, 1, NA)))
  expect_true(all(sapply(res[-3], is.vector)))
  expect_true(all(sapply(res[-3], function(x) !is.null(names(x)))))
})
