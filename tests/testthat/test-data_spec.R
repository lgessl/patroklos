test_that("prepend_to_directory() works", {

  model_1 <- ModelSpec("model", "mock1", zeroSum::zeroSum, split_index = 1, time_cutoffs = 2)
  model_2 <- model_1
  model_2$directory <- "mock2"
  model_list <- list(model_1, model_2)

  expect_silent(new_msl <- prepend_to_directory(model_list, "prepend_me"))
  expect_equal(new_msl[[1]]$directory, "prepend_me/mock1")
  expect_equal(new_msl[[2]]$directory, "prepend_me/mock2")
})
