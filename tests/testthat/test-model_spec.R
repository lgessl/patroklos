test_that("at_time_cutoff() works", {

  model <- ModelSpec(
    name = "cox",
    directory = "models/cox",
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
  )

  model_cutoff <- at_time_cutoff(model, 1.5)
  expect_equal(model_cutoff$name, "cox@1.5")
  expect_equal(model_cutoff$time_cutoffs, 1.5)
  expect_equal(model_cutoff$directory, "models/cox/1-5")
  
  expect_error(at_time_cutoff(model, 5))
})
