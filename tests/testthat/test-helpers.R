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

test_that("intersect_by_names() works", {

  a <- 1:3
  names(a) <- 1:3
  b <- 1:3
  names(b) <- 2:4

  # Vectors
  # without NA
  a_b <- intersect_by_names(a, b)
  expect_equal(names(a_b[[1]]), names(a_b[[2]]))
  # with NA
  a[2] <- NA
  a_b <- intersect_by_names(a, b, rm_na = TRUE)
  expect_equal(names(a_b[[1]]), names(a_b[[2]]))
  expect_length(a_b[[1]], 1)

  # Matrices
  # without NA
  a <- as.matrix(a)
  b <- as.matrix(b)
  a_b <- intersect_by_names(a, b)
  expect_equal(rownames(a_b[[1]]), c("2", "3"))
  expect_equal(rownames(a_b[[2]]), c("2", "3"))
  # with NA
  a[2, 1] <- NA
  a_b <- intersect_by_names(a, b, rm_na = TRUE)
  expect_equal(rownames(a_b[[2]]), c("3"))
  expect_equal(rownames(a_b[[1]]), rownames(a_b[[2]]))

  # Errors
  # input types
  expect_error(intersect_by_names(a, b[1, ]), "both matrices or")
  # matrices
  rownames(a) <- 5:7
  expect_error(intersect_by_names(a, b), "no common row names")
  rownames(a) <- NULL
  expect_error(intersect_by_names(a, b), "a has no row names")
  # vectors
  a <- a[, 1]
  b <- b[, 1]
  expect_error(intersect_by_names(a, b), "no names")
  names(a) <- 5:7
  expect_error(intersect_by_names(a, b), "no common names")
})

test_that("mirror_directory() works", {

  fp <- file.path("path", "to", "some_file")
  expect_equal(
    mirror_directory(fp, c("path", "new_path")),
    file.path("new_path", "to", "some_file")
  )
  expect_equal(
    mirror_directory(fp, c("to","from")),
    file.path("path", "from", "some_file")
  )
  fp <- stringr::str_replace(fp, "some_file", "to")
  expect_error(
    mirror_directory(fp, c("to","from")),
    "exactly one element"
  )
})
