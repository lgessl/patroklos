test_that("tibble_to_binary works correctly", {
  # Create a tibble for testing
  test_tibble <- tibble::tibble(
    column1 = factor(c("a", "b", "c")),
    column2 = factor(c("d", "e", "f"))
  )
  result <- tibble_to_binary(test_tibble)

  # Check the result
  expect_equal(
    dim(result), 
    c(3, 4)
    ) # 3 rows, 4 columns (2 original columns, each with 3-1 levels)
  expect_equal(colnames(result), c("column1_b", "column1_c", "column2_e", "column2_f"))
  expect_true(all(result %in% c(0., 1.)))
  expect_error(tibble_to_binary(test_tibble[["column1"]]))
})
