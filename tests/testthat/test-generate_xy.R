test_that("generate_predictor works", {

  expr_mat <- matrix(seq(1, 12, by = 1.), nrow = 3, ncol = 4)
  rownames(expr_mat) <- stringr::str_c("patient_", 1:3)
  colnames(expr_mat) <- stringr::str_c("gene_", 1:4)
  pheno_df <- tibble::tibble(
    patient_id = rownames(expr_mat),
    continuous_var = c(1, 2, 3), # +1 column
    discrete_var = c("A", "B", "A") # +1 column
  )
  include_from_continuous_pheno <- "continuous_var"
  include_from_discrete_pheno <- "discrete_var"
  
  # Test case 1: include both continuous and discrete pheno variables
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    patient_id_col = "patient_id",
    include_from_continuous_pheno,
    include_from_discrete_pheno
  )
  expect_identical(dim(x), c(3L, 6L))
  expect_identical(
    colnames(x)[5:6], 
    c("continuous_var", "discrete_var_B")
  )
  expect_identical(colnames(x)[1:4], colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Test case 2: include only continuous pheno variables
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    patient_id_col = "patient_id",
    include_from_continuous_pheno,
    NULL
  )
  expect_identical(dim(x), c(3L, 5L))
  expect_identical(
    colnames(x)[5], 
    "continuous_var"
  )
  expect_identical(colnames(x)[1:4], colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Test case 3: include only discrete pheno variables
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    patient_id_col = "patient_id",
    NULL,
    include_from_discrete_pheno
  )
  expect_identical(dim(x), c(3L, 5L))
  expect_identical(
    colnames(x)[5], 
    "discrete_var_B"
  )
  expect_identical(colnames(x)[1:4], colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")

  # Test case 4: include no pheno variables
  x <- generate_predictor(
    expr_mat = expr_mat,
    pheno_tbl = pheno_df,
    patient_id_col = "patient_id",
    NULL,
    NULL
  )
  expect_identical(dim(x), c(3L, 4L))
  expect_identical(colnames(x), colnames(expr_mat))
  expect_identical(rownames(x), rownames(expr_mat))
  expect_type(x, "double")
})

# Define test cases
test_that("generate_response works", {
  # Create a sample pheno_tbl
  pheno_tbl <- tibble::tibble(
    "pfs" = c(1.5, 2.5, 3.0, 4.0),
    "prog" = c(0, 1, 0, 0),
    "use_column" = c("A", "B", "C", "D"),
    "patient" = 1:4
  )
  
  # Test case 1: model == "lasso_zerosum"
  model <- "lasso_zerosum"
  pfs_leq <- 1.86
  y <- generate_response(
    pheno_tbl = pheno_tbl,
    model = model,
    patient_id_col = "patient",
    progression_col = "prog",
    pfs_col = "pfs",
    pfs_leq = pfs_leq 
  )
  y_expected <- matrix(c(NA, 0, 0, 0), ncol = 1)
  rownames(y_expected) <- pheno_tbl[["patient"]]
  colnames(y_expected) <- "pfs_leq_1.9"
  expect_equal(y, y_expected)
  expect_type(y, "double")
  
  # Test case 2: model == "cox_lasso_zerosum"
  model <- "cox_lasso_zerosum"
  y <- generate_response(
    pheno_tbl = pheno_tbl, 
    model = model,
    patient_id_col = "patient", 
    pfs_col = "pfs",
    progression_col = "prog",
    pfs_leq = pfs_leq
  )
  expect_equal(rownames(y), as.character(pheno_tbl[["patient"]]))
  expect_equal(dim(y), c(4L, 2L))
  expect_equal(colnames(y), c("pfs", "prog"))
  expect_type(y, "double")
})
