test_that("prepare_x() works", {

  set.seed(349)

  data <- generate_mock_data(n_samples = 50, n_genes = 4, n_na_in_pheno = 10)
  model <- Model$new(
    name = "zerosum",
    directory = "some_dir",
    fitter = ptk_zerosum,
    split_index = 1,
    time_cutoffs = 2.,
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    li_var_suffix = "APP"
  )
  pheno_tbl <- data$pheno_tbl
  
  # Include both continuous and discrete pheno variables
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

  # Include no pheno variables
  model$split_index <- 3
  model$include_from_continuous_pheno <- NULL
  model$include_from_discrete_pheno <- NULL
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_identical(colnames(x), colnames(data$expr_mat))
  rnames <- pheno_tbl[["patient_id"]][pheno_tbl[["split_3"]] == data$cohort]
  expect_identical(rownames(x), rnames)
  expect_type(x, "double")

  # Pheno variables without expression data
  model$include_expr <- FALSE
  model$include_from_continuous_pheno <- "continuous_var"
  model$include_from_discrete_pheno <- "discrete_var"
  x <- prepare_x(data = data, model = model, quiet = TRUE)
  expect_equal(ncol(x), 3)

  # Error: include categorical pheno variable as continuous
  data$pheno_tbl[["discrete_var"]] <- as.character(data$pheno_tbl[["discrete_var"]])
  model$include_from_continuous_pheno <- "discrete_var"
  expect_error(prepare_x(data = data, model = model, quiet = TRUE), "must be numeric")
})

test_that("prepare_y() works", {

  pheno_tbl <- tibble::tibble(
    "pfs" = c(1.5, 2.5, 3.0, 4.0),
    "prog" = c(0, 1, 0, 0),
    "use_column" = c("A", "B", "C", "D"),
    "patient" = 1:4
  )
  data <- Data$new(
    name = "Mock et al. (2023)",
    directory = "some_dir",
    train_prop = 0.8,
    patient_id_col = "patient",
    pivot_time_cutoff = 2,
    time_to_event_col = "pfs",
    event_col = "prog"
  )
  data$pheno_tbl <- pheno_tbl
  model <- Model$new(
    name = "zerosum",
    directory = "some_dir",
    fitter = ptk_zerosum,
    split_index = 1,
    time_cutoffs = 1.9,
    response_type = "binary"
  )

  # Binary response
  y <- prepare_y(data = data, model = model)
  y_expected <- matrix(c(NA, 0, 0, 0), ncol = 1)
  rownames(y_expected) <- pheno_tbl[["patient"]]
  colnames(y_expected) <- "time_cutoff_1.9"
  expect_equal(y, y_expected)
  expect_type(y, "double")

  model$time_cutoffs <- 3.5
  y_expected[, 1] <- c(NA, 1, NA, 0)
  colnames(y_expected) <- "time_cutoff_3.5"
  y <- prepare_y(data = data, model = model)
  expect_equal(y, y_expected)
  
  # survival_censored response
  model$response_type <- "survival_censored"
  model$time_cutoffs <- Inf
  y <- prepare_y(data = data, model = model)
  y_expected <- pheno_tbl[, c("pfs", "prog")] |> as.matrix()
  rownames(y_expected) <- pheno_tbl[["patient"]]
  colnames(y_expected) <- model$response_colnames
  expect_equal(y, y_expected)

  model$time_cutoffs <- 2.7
  y <- prepare_y(data = data, model = model)
  y_expected[c(3, 4), 1] <- model$time_cutoffs
  y_expected[c(3, 4), 2] <- 0
  expect_equal(y, y_expected)
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
    combined_feature_positive_ratio = 0.4)
  expect_equal(colnames(x_wide), c("var_1", "var_1&var_2", "var_2", "var_2&var_3", 
    "var_3"))
  exp <- c(0.7, 0.6, 0)
  names(exp) <- rownames(x)
  expect_equal(x_wide[, "var_2&var_3"], exp)
  expect_equal(rownames(x), rownames(x_wide))

  x[1, 1] <- NA
  x_wide <- combine_features(x, combine_n_max_features = 3, 
    combined_feature_positive_ratio = 0)
  expect_equal(colnames(x_wide), c("var_1", "var_1&var_2", "var_1&var_3",
    "var_1&var_2&var_3", "var_2", "var_2&var_3", "var_3"))
  exp <- c(NA, 0, 0)
  names(exp) <- rownames(x)
  expect_equal(x_wide[, "var_1&var_2&var_3"], exp)

  x <- as.matrix(c(1, 0, 1))
  colnames(x) <- "var_1"
  x_wide <- combine_features(x, combine_n_max_features = 1, 
    combined_feature_positive_ratio = 0.4)
  expect_equal(x, x_wide)

  x <- matrix(c(c(1,0,NA), c(0,1,1)), ncol = 2)
  colnames(x) <- c("var_1", "var_2")
  x_wide <- combine_features(x, combine_n_max_features = 1, 
    combined_feature_positive_ratio = 0.4)
  expect_equal(x, x_wide)
})
