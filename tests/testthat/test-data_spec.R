test_that("prepend_to_directory() works", {

  model_spec_list = list(
    ModelSpec("mock1", fitter = zeroSum::zeroSum, directory = "mock1"),
    ModelSpec("mock2", fitter = zeroSum::zeroSum, directory = "mock2")
  )
  expect_silent(new_msl <- prepend_to_directory(model_spec_list, "prepend_me"))
  expect_equal(new_msl[[1]]$directory, "prepend_me/mock1")
  expect_equal(new_msl[[2]]$directory, "prepend_me/mock2")
})
