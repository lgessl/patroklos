test_that("Model$new() works", {
  
    model <- Model$new(
      name = "cox",
      fitter = ptk_zerosum,
      directory = "models/cox",
      split_index = 1:2,
      time_cutoffs = c(1.5, 2)
    )
    expect_equal(model$name, "cox")
    expect_equal(model$fitter, ptk_zerosum)
    expect_equal(model$directory, "models/cox")
    expect_equal(model$split_index, 1:2)
    expect_equal(model$time_cutoffs, c(1.5, 2))
  
})

test_that("Model$at_time_cutoff() works", {

  model <- Model$new(
    name = "cox",
    fitter = ptk_zerosum,
    directory = "models/cox",
    split_index = 1:2,
    time_cutoffs = c(1.5, 2)
  )

  model_cutoff <- model$at_time_cutoff(1.5)
  expect_equal(model_cutoff$name, "cox@1.5")
  expect_equal(model_cutoff$time_cutoffs, 1.5)
  expect_equal(model_cutoff$directory, "models/cox/1-5")
  
  expect_error(model$at_time_cutoff(5))
})

test_that("prepend_to_directory() works", {

  model_1 <- Model$new(
    name = "model", 
    directory = "mock1", 
    fitter = ptk_zerosum, 
    split_index = 1, 
    time_cutoffs = 2
  )
  model_2 <- model_1$clone()
  model_2$directory <- "mock2"
  models <- list(model_1, model_2)

  new_models <- prepend_to_directory(models, "prepend_me")
  expect_equal(new_models[[1]]$directory, "prepend_me/mock1")
  expect_equal(new_models[[2]]$directory, "prepend_me/mock2")
  expect_equal(model_1$directory, "mock1")
  expect_equal(model_2$directory, "mock2")
})
