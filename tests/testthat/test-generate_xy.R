test_that("generate_predictor() works", {

  expr_mat <- matrix(seq(1, 12, by = 1.), nrow = 3, ncol = 4)
  rownames(expr_mat) <- stringr::str_c("patient_", 1:3)
  colnames(expr_mat) <- stringr::str_c("gene_", 1:4)
  pheno_df <- tibble::tibble(
    patient_id = rownames(expr_mat),
    continuous_var = c(1, 2, 3), # +1 column
    discrete_var = c("A", "B", "A") # +1 column
  )
  data <- Data$new(
    name = "Mock et al. (2023)", 
    directory =  "some_dir",
    train_prop = 0.8
  )
  model <- Model$new(
    name = "zerosum",
    directory = "some_dir",
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2.,
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    append_to_includes = "APP"
  )
  
  # Test case 1: include both continuous and discrete pheno variables
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    data = data,
    model = model
  )
  expect_identical(
    colnames(x)[5:6], 
    c("continuous_varAPP", "discrete_var_BAPP")
  )
  expect_identical(colnames(x)[1:4], colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Test case 2: include only continuous pheno variables
  model$include_from_discrete_pheno <- NULL
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    data = data,
    model = model
  )
  expect_identical(
    colnames(x)[5], 
    "continuous_varAPP"
  )
  expect_identical(colnames(x)[1:4], colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Test case 3: include only discrete pheno variables
  model$include_from_continuous_pheno <- NULL
  model$include_from_discrete_pheno <- "discrete_var"
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    data = data,
    model = model
  )
  expect_identical(dim(x), c(3L, 5L))
  expect_identical(
    colnames(x)[5], 
    "discrete_var_BAPP"
  )
  expect_identical(colnames(x)[1:4], colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Test case 4: include no pheno variables
  model$include_from_continuous_pheno <- NULL
  model$include_from_discrete_pheno <- NULL
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    data = data,
    model = model
  )
  expect_identical(colnames(x), colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Error: include categorical pheno variable as continuous
  model$include_from_continuous_pheno <- "discrete_var"
  expect_error(
    generate_predictor(
      expr_mat = expr_mat,
      pheno_tbl = pheno_df,
      data = data,
      model = model
    ),
    "must be numeric"
  )
})


test_that("generate_response() works", {

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
    time_to_event_col = "pfs",
    event_col = "prog"
  )
  model <- Model$new(
    name = "zerosum",
    directory = "some_dir",
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 1.9,
    response_type = "binary"
  )

  # Test case 1: binary response
  y <- generate_response(
    pheno_tbl = pheno_tbl,
    data = data,
    model = model
  )
  y_expected <- matrix(c(NA, 0, 0, 0), ncol = 1)
  rownames(y_expected) <- pheno_tbl[["patient"]]
  colnames(y_expected) <- "time_cutoff_1.9"
  expect_equal(y, y_expected)
  expect_type(y, "double")

  model$time_cutoffs <- 3.5
  y_expected[, 1] <- c(NA, 1, NA, 0)
  colnames(y_expected) <- "time_cutoff_3.5"
  y <- generate_response(
    pheno_tbl = pheno_tbl,
    data = data,
    model = model
  )
  expect_equal(y, y_expected)
  
  # Test case 2: survival_censored response
  model$response_type <- "survival_censored"
  model$time_cutoffs <- Inf
  y <- generate_response(
    pheno_tbl = pheno_tbl, 
    data = data,
    model = model
  )
  y_expected <- pheno_tbl[, c("pfs", "prog")] |> as.matrix()
  rownames(y_expected) <- pheno_tbl[["patient"]]
  colnames(y_expected) <- model$response_colnames
  expect_equal(y, y_expected)

  model$time_cutoffs <- 2.7
  y <- generate_response(
    pheno_tbl = pheno_tbl, 
    data = data,
    model = model
  )
  y_expected[c(3, 4), 1] <- model$time_cutoffs
  y_expected[c(3, 4), 2] <- 0
  expect_equal(y, y_expected)
})
