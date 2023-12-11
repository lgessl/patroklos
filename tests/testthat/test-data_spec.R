test_that("prepend_to_save_dir() works", {

  model_spec_list = list(
    ModelSpec("mock1", fitter = zeroSum::zeroSum, save_dir = "mock1"),
    ModelSpec("mock2", fitter = zeroSum::zeroSum, save_dir = "mock2")
  )
  expect_silent(new_msl <- prepend_to_save_dir(model_spec_list, "prepend_me"))
  expect_equal(new_msl[[1]]$save_dir, "prepend_me/mock1")
  expect_equal(new_msl[[2]]$save_dir, "prepend_me/mock2")
})
