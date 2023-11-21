test_that("generate_predictor works", {

  expr_mat <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(expr_mat) <- stringr::str_c("patient_", 1:3)
  pheno_df <- tibble::tibble(
    continuous_var = c(1, 2, 3), # +1 column
    discrete_var = c("A", "B", "A") # +1 column
  )
  include_from_continuous_pheno <- "continuous_var"
  include_from_discrete_pheno <- "discrete_var"
  
  x <- generate_predictor(
    expr_mat,
    pheno_df,
    include_from_continuous_pheno,
    include_from_discrete_pheno
  )
  
  expect_identical(dim(x), c(3L, 6L))
  expect_identical(
    colnames(x)[5:6], 
    c("continuous_var", "discrete_var_B")
  )
  expect_identical(
    rownames(x), rownames(expr_mat)
  )
})

# Define test cases
test_that("generate_response works", {
  # Create a sample pheno_tbl
  pheno_tbl <- tibble::tibble(
    pfs_yrs = c(1.5, 2.5, 3.0, 4.0),
    progression = c(0, 1, 0, 0),
    use_column = c("A", "B", "C", "D"),
    "patient_id" = 1:4
  )
  
  # Test case 1: model == "lasso_zerosum"
  model <- "lasso_zerosum"
  pfs_leq <- 2.0
  expected <- rep(FALSE, 3)
  names(expected) <- c("2", "3", "4")
  output <- generate_response(pheno_tbl, model, pfs_leq)
  expect_equal(output, expected)
  
  # Test case 2: model == "cox_lasso_zerosum"
  model <- "cox_lasso_zerosum"
  expected <- pheno_tbl[, c("pfs_yrs", "progression")]
  output <- generate_response(pheno_tbl, model)
  expect_equal(output, expected)
})
