test_that("Model$new() works", {
  
    model <- Model$new(
      name = "cox",
      fitter = zeroSum::zeroSum,
      directory = "models/cox",
      split_index = 1:2,
      time_cutoffs = c(1.5, 2)
    )
    expect_equal(model$name, "cox")
    expect_equal(model$fitter, zeroSum::zeroSum)
    expect_equal(model$directory, "models/cox")
    expect_equal(model$split_index, 1:2)
    expect_equal(model$time_cutoffs, c(1.5, 2))
  
})

test_that("Model$at_time_cutoff() works", {

  model <- Model$new(
    name = "cox",
    fitter = zeroSum::zeroSum,
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
