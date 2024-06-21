test_that("Model$new() works", {
  
    model <- Model$new(
      name = "cox",
      fitter = ptk_zerosum,
      directory = "models/cox",
      split_index = 1:2,
      time_cutoffs = c(1.5, 2),
      val_error_fun = neg_roc_auc,
      continuous_output = TRUE
    )
    expect_equal(model$name, "cox")
    expect_equal(model$fitter, ptk_zerosum)
    expect_equal(model$directory, "models/cox")
    expect_equal(model$split_index, 1:2)
    expect_equal(model$time_cutoffs, c(1.5, 2))
  
})

test_that("prepend_to_directory() works", {

  model1 <- Model$new(
    name = "model1", 
    directory = "mock1", 
    fitter = ptk_zerosum, 
    split_index = 1, 
    time_cutoffs = 2,
    val_error_fun = neg_roc_auc,
    continuous_output = TRUE
  )
  model2 <- model1$clone()
  model2$name <- "model2"
  model2$directory <- "mock2"
  models <- list(model1, model2)

  prepend_to_directory(models, "prepend_me")
  expect_equal(model1$directory, "prepend_me/mock1")
  expect_equal(model2$directory, "prepend_me/mock2")

  models[[2]]$directory <- "prepend_me/mock1"
  expect_error(prepend_to_directory(models, "prepend_me"), 
    "directories are not unique")
  models[[2]]$directory <- "mock2"

  models[[2]]$name <- "model1"
  expect_error(prepend_to_directory(models, "prepend_me"), 
    "names are not unique")
})
