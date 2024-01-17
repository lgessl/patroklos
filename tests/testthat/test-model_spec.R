test_that("at_time_cutoff() works", {

  model_spec <- ModelSpec(
    name = "cox",
    directory = "models/cox",
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = c(1.5, 2),
  )

  ms_cutoff <- at_time_cutoff(model_spec, 1.5)
  expect_equal(ms_cutoff$name, "cox@1.5")
  expect_equal(ms_cutoff$time_cutoffs, 1.5)
  expect_equal(ms_cutoff$directory, "models/cox/1-5")
  
  expect_error(at_time_cutoff(model_spec, 5))
})
