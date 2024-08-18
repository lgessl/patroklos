test_that("Data$new() works", {

  data <- Data$new(
    name = "mock",
    directory = "mock",
    pivot_time_cutoff = 0.5,
    expr_file = "expr.csv",
    pheno_file = "pheno.csv",
    cohort = "train",
    patient_id_col = "patient_id",
    time_to_event_col = "time_to_event",
    event_col = "event",
    cohort_col = "split_1",
    benchmark_col = "benchmark",
    gene_id_col = "gene_id",
  )
  expect_equal(data$name, "mock")
  expect_equal(data$directory, "mock")
  expect_equal(data$pivot_time_cutoff, 0.5)
  expect_equal(data$expr_file, "expr.csv")
  expect_equal(data$pheno_file, "pheno.csv")
  expect_equal(data$cohort, "train")
  expect_equal(data$patient_id_col, "patient_id")
  expect_equal(data$time_to_event_col, "time_to_event")
  expect_equal(data$event_col, "event")
  expect_equal(data$cohort_col, "split_1")
  expect_equal(data$benchmark_col, "benchmark")
  expect_equal(data$gene_id_col, "gene_id")
  expect_equal(data$pivot_time_cutoff, 0.5)
})

test_that("Data$read() works", {

  n_samples <- 5
  n_genes <- 2

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = 3,
    to_csv = dir
  ) 
  ncol_pheno <- ncol(data$pheno_tbl)
  data$expr_mat <- NULL
  data$pheno_tbl <- NULL
  data$read()

  expect_true(is.numeric(data$expr_mat))
  expect_true(is.matrix(data$expr_mat))
  expect_s3_class(data$pheno_tbl, "tbl_df")
  expect_equal(colnames(data$pheno_tbl)[1], data$patient_id_col)
  expect_equal(colnames(data$expr_mat), paste0("gene_", 1:n_genes))
  expect_equal(rownames(data$expr_mat), paste0("sample_", 1:n_samples))
  expect_equal(data$pheno_tbl[[data$patient_id_col]], paste0("sample_", 1:n_samples))
  expect_equal(dim(data$expr_mat), c(n_samples, n_genes))
  expect_equal(dim(data$pheno_tbl), c(n_samples, ncol_pheno))

  # errors
  data$patient_id_col <- "missing_patient_id_col"
  expect_error(data$read(), "does not exist in")
  data$gene_id_col <- "missing_gene_id"
  expect_error(data$read(), "First column of expression")
})

test_that("Data$survival_quantiles()", {
  
  set.seed(321)

  n_samples <- 20

  data <- generate_mock_data(n_samples = n_samples, n_genes = 1)
  
  tbl <- data$survival_quantiles()
  expect_true(nrow(tbl) <= n_samples)
  q <- tbl[["quantile"]]
  v <- tbl[[data$time_to_event_col]]
  expect_true(all(q >= 0 & q <= 1))
  expect_true(all(v >= 0))
  expect_equal(order(q), seq_along(q))
  expect_equal(order(v), seq_along(v))
})

test_that("Data$split() works", {

  set.seed(234)

  n_samples <- 100
  n_genes <- 1

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = 0,
  )
  data$directory <- withr::local_tempdir()
  data$pheno_tbl[["split_1"]] <- NULL

  pt_old <- data$pheno_tbl
  data$split(train_prop = 0.8, save = TRUE, keep_risk = TRUE)
  expect_equal(nrow(pt_old), nrow(data$pheno_tbl))
  expect_equal(ncol(pt_old)+1, ncol(data$pheno_tbl))
  expect_false(is.null(data$pheno_tbl[["split_1"]]))
  expect_true(file.exists(file.path(data$directory, "cohort.rds")))

  expect_message(data$split(train_prop = 0.8, save = TRUE), 
    "Found a cohort file")

  cohort <- readRDS(file.path(data$directory, "cohort.rds"))
  saveRDS(cohort[-1], file.path(data$directory, "cohort.rds"))
  expect_error(data$split(train_prop = 0.7, save = TRUE, quiet = TRUE), 
    "The patient IDs in the cohort file do not match")
  unlink(file.path(data$directory, "cohort.rds"))

  data$split(train_prop = 0.7, save = FALSE, keep_risk = FALSE)
  expect_false(file.exists(file.path(data$directory, "cohort.rds")))
})
