test_that("check_fitter() works correctly", {

  valid_fitter <- zeroSum::zeroSum
  invalid_fitter1 <- function(y, ...) {}
  invalid_fitter2 <- function(x) {}
  
  expect_silent(check_fitter(valid_fitter))

  optional_fitter_args <- list("alpha" = 0.5, nFold = 3)
  expect_silent(check_fitter(valid_fitter, optional_fitter_args))

  wrong_fitter_args <- list("alpha" = 0.5, nFold = 3, "wrong_arg" = 1)
  expect_silent(check_fitter(valid_fitter, wrong_fitter_args))

  expect_error(check_fitter(invalid_fitter1))
  expect_error(check_fitter(invalid_fitter2))

})