test_that("check_fitter() works correctly", {

  valid_fitter1 <- ptk_zerosum
  valid_fitter2 <- function(x, y_bin, y_cox, beta) {}
  invalid_fitter <- function(y, ...) {}
  
  check_fitter(valid_fitter1)

  hyperparams <- list(beta = 0.5)
  check_fitter(valid_fitter2, hyperparams)

  wrong_fitter_args <- list(alpha = 0.5)
  expect_error(check_fitter(valid_fitter, wrong_fitter_args))

  expect_error(check_fitter(invalid_fitter))
})