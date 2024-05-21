test_that("Data$new() works", {

  data <- Data$new(
    name = "mock",
    directory = "mock",
    train_prop = 0.7,
    pivot_time_cutoff = 0.5,
    expr_file = "expr.csv",
    pheno_file = "pheno.csv",
    cohort = "train",
    patient_id_col = "patient_id",
    time_to_event_col = "time_to_event",
    event_col = "event",
    benchmark_col = "benchmark",
    gene_id_col = "gene_id",
    split_col_prefix = "split"
  )
  expect_equal(data$name, "mock")
  expect_equal(data$directory, "mock")
  expect_equal(data$train_prop, 0.7)
  expect_equal(data$pivot_time_cutoff, 0.5)
  expect_equal(data$expr_file, "expr.csv")
  expect_equal(data$pheno_file, "pheno.csv")
  expect_equal(data$cohort, "train")
  expect_equal(data$patient_id_col, "patient_id")
  expect_equal(data$time_to_event_col, "time_to_event")
  expect_equal(data$event_col, "event")
  expect_equal(data$benchmark_col, "benchmark")
  expect_equal(data$gene_id_col, "gene_id")
  expect_equal(data$split_col_prefix, "split")
  expect_equal(data$pivot_time_cutoff, 0.5)
})

test_that("Data$read() works", {

  n_samples <- 5
  n_genes <- 2

  dir <- withr::local_tempdir()
  generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = 3,
    to_csv = dir,
    split_index = 1:3
  ) 
  data <- Data$new(
    name = "mock",
    directory = dir,
    train_prop = 0.7,
    pivot_time_cutoff = 0.5,
    expr_file = "expr.csv",
    pheno_file = "pheno.csv",
    cohort = "train",
    patient_id_col = "patient_id",
    time_to_event_col = "time_to_event",
    event_col = "event",
    benchmark_col = "benchmark",
    gene_id_col = "gene_id",
    split_col_prefix = "split"
  )
  data$read()

  expect_true(is.numeric(data$expr_mat))
  expect_true(is.matrix(data$expr_mat))
  expect_s3_class(data$pheno_tbl, "tbl_df")
  expect_equal(colnames(data$pheno_tbl)[1], data$patient_id_col)
  expect_equal(colnames(data$expr_mat), paste0("gene_", 1:n_genes))
  expect_equal(rownames(data$expr_mat), paste0("sample_", 1:n_samples))
  expect_equal(data$pheno_tbl[[data$patient_id_col]], paste0("sample_", 1:n_samples))

  # errors
  # directory does not exist
  data$patient_id_col <- "missing_patient_id_col"
  expect_error(data$read(), "no column named")
  data$gene_id_col <- "missing_gene_id"
  expect_error(data$read(), "no column named")
})

test_that("Data$survival_quantiles()", {
  
  set.seed(321)

  n_samples <- 20

  data <- generate_mock_data(n_samples = n_samples, n_genes = 1, split_index = 1)
  
  q <- data$survival_quantiles()
  expect_true(length(q) <= n_samples)
  expect_true(all(names(q) >= 0 & names(q) <= 1))
  expect_true(all(q >= 0))
  expect_equal(order(q), seq_along(q))
})
