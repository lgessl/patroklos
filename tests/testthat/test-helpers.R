test_that("dichotomize_tibble works correctly", {
  # Create a tibble for testing
  test_tibble <- tibble::tibble(
    column1 = factor(c("a", "b", "b", "c")),
    column2 = factor(c("d", "e", "f", "d"))
  )
  level_list <- lapply(test_tibble, function(c) levels(as.factor(c)))
  result <- dichotomize_tibble(test_tibble, level_list = level_list)

  # Check the result
  expect_equal(
    dim(result), 
    c(4, 4)
    ) # 3 rows, 4 columns (2 original columns, each with 3-1 levels)
  expect_equal(colnames(result), c("column1_b", "column1_c", "column2_e", "column2_f"))
  expect_true(all(result %in% c(0., 1.)))
  expect_error(dichotomize_tibble(test_tibble[["column1"]]))

  test_tibble[1, 2] <- NA
  result <- dichotomize_tibble(test_tibble, level_list = level_list)
  expect_true(all(is.na(result[1, 3:4])))
})

test_that("intersect_by_names() works", {

  a <- 1:4
  names(a) <- 1:4
  b <- 1:4
  names(b) <- 2:5

  # Vectors
  # without NA
  a_b <- intersect_by_names(a, b)
  expect_equal(names(a_b[[1]]), names(a_b[[2]]))
  # with NA
  a[2] <- NA
  a_b <- intersect_by_names(a, b, rm_na = c(TRUE, TRUE))
  expect_equal(names(a_b[[1]]), names(a_b[[2]]))
  expect_length(a_b[[1]], 2)

  # Matrices
  # without NA
  a <- as.matrix(a)
  b <- as.matrix(b)
  b[2,1] <- NA
  a_b <- intersect_by_names(a, b, rm_na = c(TRUE, FALSE))
  expect_equal(rownames(a_b[[1]]), c("3", "4"))
  expect_equal(rownames(a_b[[2]]), c("3", "4"))
  # with NA
  a[2, 1] <- NA
  a_b <- intersect_by_names(a, b, rm_na = c(TRUE, TRUE))
  expect_equal(rownames(a_b[[2]]), c("4"))
  expect_equal(rownames(a_b[[1]]), rownames(a_b[[2]]))

  # Errors
  # input types
  expect_error(intersect_by_names(a, b[1, ]), "both matrices or")
  # matrices
  rownames(a) <- 10:13
  expect_error(intersect_by_names(a, b), "no common row names")
  rownames(a) <- NULL
  expect_error(intersect_by_names(a, b), "a has no row names")
  # vectors
  a <- a[, 1]
  b <- b[, 1]
  expect_error(intersect_by_names(a, b), "no names")
  names(a) <- 5:8
  names(b) <- 1:4
  expect_error(intersect_by_names(a, b), "no common names")
})

test_that("create_data_partition() works", {

  n <- 50
  n_levels <- 3

  outcome <- sample(1:n_levels, n, replace = TRUE)
  outcome <- as.factor(outcome)

  for(p in c(.2, .8)){
    selected <- create_data_partition(outcome, p = p)
    expect_true(is.numeric(selected))
    expect_true(all(
      length(selected) >= p*n - n_levels/2,
      length(selected) <= p*n + n_levels/2
    ))
  }
})