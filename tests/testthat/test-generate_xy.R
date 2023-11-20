test_that("generate_predictor works", {

  expr_df <- matrix(1:12, nrow = 4)
  pheno_df <- tibble::tibble(
    continuous_var = c(1, 2, 3),
    discrete_var = c("A", "B", "A")
  )
  include_from_continuous_pheno <- "continuous_var"
  include_from_discrete_pheno <- "discrete_var"
  
  result <- generate_predictor(
    expr_df,
    pheno_df,
    include_from_continuous_pheno,
    include_from_discrete_pheno
  )
  
  expect_identical(dim(result), c(3L, 6L))
  expect_identical(
    colnames(result)[5:6], 
    c("continuous_var", "discrete_var_B")
  )
})

# Define test cases
test_that("generate_response works", {
  # Create a sample pheno_tbl
  pheno_tbl <- tibble::tibble(
    pfs_yrs = c(1.5, 2.5, 3.0, 4.0),
    progression = c(0, 1, 0, 0),
    use_column = c("A", "B", "C", "D")
  )
  rownames(pheno_tbl) <- 1:4
  
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
