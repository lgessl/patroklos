test_that("prepare() works", {

  set.seed(343)

  data <- generate_mock_data(n_samples = 50, n_genes = 4, n_na_in_pheno = 10)
  model <- Model$new(
    name = "cox",
    fitter = ptk_zerosum,
    directory = "some_dir",
    time_cutoffs = Inf,
    val_error_fun = neg_roc_auc,
    hyperparams = list(family = "cox", zeroSum = FALSE),
    include_from_continuous_pheno = c("continuous_var"),
    include_from_discrete_pheno = c("discrete_var", "abc_gcb", "ipi_age"),
    combine_n_max_categorical_features = 1:2,
    combined_feature_min_positive_ratio = 0.04
  )
  xy <- data$prepare(model, quiet = TRUE)
  expect_true(any(stringr::str_detect(colnames(xy[[1]]), "&")))
  expect_equal(ncol(xy[[2]]), 2)
  expect_equal(length(unique(xy[[2]][, 1])), nrow(xy[[2]]))
  expect_true(all(xy[[2]][, 2] %in% c(0, 1)))
})

test_that("prepare_x() works", {

  set.seed(349)

  data <- generate_mock_data(n_samples = 50, n_genes = 4, n_na_in_pheno = 10)
  model <- Model$new(
    name = "zerosum",
    directory = "some_dir",
    fitter = ptk_zerosum,
    time_cutoffs = 2.,
    val_error_fun = neg_roc_auc,
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    li_var_suffix = "APP"
  )
  pheno_tbl <- data$pheno_tbl
  
  # Include both continuous and discrete pheno variables
  data$cohort <- "train"
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_identical(
    colnames(x)[5:7], 
    c("continuous_varAPP", "discrete_var_2APP", "discrete_var_3APP")
  )
  expect_identical(colnames(x)[1:4], colnames(data$expr_mat))
  rnames <- pheno_tbl[["patient_id"]][pheno_tbl[["split_1"]] == data$cohort]
  expect_true(all(rownames(x) %in% rnames))
  expect_type(x, "double")

  # Less samples, mean impute
  data$cohort <- "test"
  data$expr_mat <- data$expr_mat[5:35, ]
  data$pheno_tbl <- data$pheno_tbl[5:35, ]
  data$imputer <- mean_impute
  x_small <- prepare_x(data = data, model = model, quiet = TRUE)
  pheno_tbl <- data$pheno_tbl
  rnames <- pheno_tbl[["patient_id"]][pheno_tbl[["split_1"]] == data$cohort]
  expect_identical(colnames(x), colnames(x_small))
  expect_equal(rownames(x_small), rnames)

  # Include only continuous pheno variables
  model$include_from_discrete_pheno <- NULL
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_identical(colnames(x)[5], "continuous_varAPP")
  expect_identical(colnames(x)[1:4], colnames(data$expr_mat))
  rnames <- pheno_tbl[["patient_id"]][pheno_tbl[["split_1"]] == data$cohort]
  expect_identical(rownames(x), rnames)
  expect_type(x, "double")

  # Include only discrete pheno variables
  model$include_from_continuous_pheno <- NULL
  model$include_from_discrete_pheno <- "discrete_var"
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_equal(ncol(x), ncol(data$expr_mat)+2)
  expect_identical(colnames(x)[5], "discrete_var_2APP")
  expect_identical(colnames(x)[1:4], colnames(data$expr_mat))
  rnames <- pheno_tbl[["patient_id"]][pheno_tbl[["split_1"]] == data$cohort]
  expect_identical(rownames(x), rnames)
  expect_type(x, "double")

  model$combine_n_max_categorical_features <- 1:2
  model$combined_feature_min_positive_ratio <- 0
  model$include_from_discrete_pheno <- c("discrete_var", "abc_gcb")

  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_equal(ncol(x), 4+4+4)

  # Include no pheno variables
  model$include_from_discrete_pheno <- NULL
  model$include_from_continuous_pheno <- NULL
  model$include_from_discrete_pheno <- NULL
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_identical(colnames(x), colnames(data$expr_mat))
  rnames <- pheno_tbl[["patient_id"]][pheno_tbl[["split_1"]] == data$cohort]
  expect_identical(rownames(x), rnames)
  expect_type(x, "double")

  # Pheno variables without expression data
  model$include_expr <- FALSE
  model$include_from_continuous_pheno <- "continuous_var"
  model$include_from_discrete_pheno <- "discrete_var"
  model$enable_imputation <- FALSE
  data$pheno_tbl[1, "continuous_var"] <- NA
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_equal(ncol(x), 3)
  expect_true(is.na(x[1, 1]))
  model$enable_imputation <- TRUE

  # Error: include categorical pheno variable as continuous
  data$pheno_tbl[["discrete_var"]] <- as.character(data$pheno_tbl[["discrete_var"]])
  model$include_from_continuous_pheno <- "discrete_var"
  expect_error(prepare_x(data = data, model = model, quiet = TRUE), "must be numeric")

  model$include_from_continuous_pheno <- NULL
  data$cohort <- "test|train"
  expect_error(prepare_x(data = data, model = model, quiet = TRUE), 
    "All patients are in the same cohort")
})

test_that("binarize_y() works", {

  y_cox <- matrix(c(c(1.5, 2.5, 1, 3, 2.1), c(0, 1, 1, 0, 1)), ncol = 2)
  rownames(y_cox) <- paste0("sample_", 1:nrow(y_cox))
  y_bin <- binarize_y(y_cox, time_cutoff = 2, pivot_time_cutoff = 2)
  y_exp <- as.matrix(c(NA, 0, 1, 0, 0))
  rownames(y_exp) <- rownames(y_cox)
  expect_equal(y_bin, y_exp)

  y_bin <- binarize_y(y_cox, time_cutoff = Inf, pivot_time_cutoff = 2)
  expect_equal(y_bin, y_exp)

  y_bin <- binarize_y(y_cox, time_cutoff = 2.2, pivot_time_cutoff = 2)
  y_exp[, 1] <- c(NA, 0, 1, 0, 1)
  expect_equal(y_bin, y_exp)
})

test_that("censor_y() works", {

  y_cox <- matrix(c(c(1.5, 2.5, 1, 3, 2.1), c(0, 1, 1, 0, 1)), ncol = 2)
  rownames(y_cox) <- paste0("sample_", 1:nrow(y_cox))
  y_cox_cens <- censor_y(y_cox, time_cutoff = 2)
  y_exp <- matrix(c(c(1.5, 2, 1, 2, 2), c(0, 0, 1, 0, 0)), ncol = 2)
  rownames(y_exp) <- rownames(y_cox)
  expect_equal(y_cox_cens, y_exp)

  y_cox_cens <- censor_y(y_cox, time_cutoff = 2.2)
  y_exp[, 1] <- c(1.5, 2.2, 1, 2.2, 2.1)
  y_exp[, 2] <- c(0, 0, 1, 0, 1)
  expect_equal(y_cox_cens, y_exp)
})

test_that("mean_impute() works", {

  set.seed(123)

  n_continuous <- 10
  n_dichotomous <- 5
  n_samples <- 10
  n_na <- 10

  x <- matrix(rnorm(n_continuous*n_samples), nrow = n_samples)
  x <- cbind(x, matrix(sample(0:1, n_dichotomous*n_samples, replace = TRUE), 
    nrow = n_samples))
  for (i in seq(n_na)) 
    x[sample(n_samples, 1), sample(ncol(x), 1)] <- NA
  rownames(x) <- paste0("sample_", 1:n_samples)
  colnames(x) <- paste0("var_", seq(ncol(x)))

  x_imp <- mean_impute(x)
  expect_true(all(!is.na(x_imp)))
  expect_equal(dim(x), dim(x_imp))
  expect_equal(colnames(x), colnames(x_imp))
  expect_equal(rownames(x), rownames(x_imp))
  expect_equal(x_imp[!is.na(x)], x[!is.na(x)])
  expect_true(is.matrix(x_imp))
})

test_that("combine_features() works", {

  set.seed(134)

  x <- matrix(c(c(1,0,1), c(0.7, 1, 1), c(1,0.6, 0)), ncol = 3)
  colnames(x) <- c("var_1", "var_2", "var_3")
  rownames(x) <- paste0("sample_", seq(nrow(x)))
  x_wide <- combine_features(x, combine_n_max_features = 3, 
    combined_feature_positive_ratio = 0.4, original_cnames = c("var_1", "var_2", "var_3"))
  expect_equal(colnames(x_wide), c("var_1", "var_1&var_2", "var_2", "var_2&var_3", 
    "var_3"))
  exp <- c(0.7, 0.6, 0)
  names(exp) <- rownames(x)
  expect_equal(x_wide[, "var_2&var_3"], exp)
  expect_equal(rownames(x), rownames(x_wide))

  x[1, 1] <- NA
  x_wide <- combine_features(x, combine_n_max_features = 3, 
    combined_feature_positive_ratio = 0, original_cnames = c("var_1", "var_2", "var_3"))
  expect_equal(colnames(x_wide), c("var_1", "var_1&var_2", "var_1&var_3",
    "var_1&var_2&var_3", "var_2", "var_2&var_3", "var_3"))
  exp <- c(NA, 0, 0)
  names(exp) <- rownames(x)
  expect_equal(x_wide[, "var_1&var_2&var_3"], exp)

  x <- as.matrix(c(1, 0, 1))
  colnames(x) <- "var_1"
  x_wide <- combine_features(x, combine_n_max_features = 1, 
    combined_feature_positive_ratio = 0.4, original_cnames = "var_1")
  expect_equal(x, x_wide)

  x <- matrix(c(c(1,0,NA), c(0,1,1)), ncol = 2)
  colnames(x) <- c("var_1", "var_2")
  x_wide <- combine_features(x, combine_n_max_features = 1, 
    combined_feature_positive_ratio = 0.4, original_cnames = c("var_1", "var_2"))
  expect_equal(x, x_wide)

  x <- matrix(c(c(0.5,0,1), c(0.5,1,0), c(1,0,1)), ncol = 3)
  colnames(x) <- c("a_1", "a_2", "b_1")
  x_wide <- combine_features(x, combine_n_max_features = 3, 
    combined_feature_positive_ratio = 0, original_cnames = c("a", "b"))
  expect_equal(colnames(x_wide), c("a_1", "a_1&b_1", "a_2", "a_2&b_1", "b_1"))
})
