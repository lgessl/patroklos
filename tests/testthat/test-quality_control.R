test_that("qc_preprocess() works", {

  set.seed(234)

  n_samples <- 5
  n_genes <- 1

  dir <- withr::local_tempdir()
  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = 0,
    to_csv = dir
  )
  expr_tbl <- data$expr_mat |> t() |> tibble::as_tibble(rownames = "gene_id")
  data$expr_mat <- NULL

  data$qc_preprocess(expr_tbl = expr_tbl)
  data$pheno_tbl[["patient_id"]] <- 1:nrow(data$pheno_tbl)
  names(expr_tbl) <- c("gene_id", 1:nrow(data$pheno_tbl))
  data$qc_preprocess(expr_tbl = expr_tbl)
  expect_true(is.character(data$pheno_tbl[["patient_id"]]))

  data$pheno_tbl[["progression"]][1] <- 2
  expect_error(data$qc_preprocess(expr_tbl = expr_tbl), regexp = "either 1")
  data$pheno_tbl[["progression"]][1] <- 0

  data$pheno_tbl[["patient_id"]][1] <- "2"
  expect_error(data$qc_preprocess(expr_tbl = expr_tbl), regexp = "Patient ids")
  data$pheno_tbl[["patient_id"]][1] <- "1"

  data$pheno_tbl[["pfs_years"]][1] <- NA
  expect_warning(data$qc_preprocess(expr_tbl = expr_tbl), regexp = "missing values")
  data$pheno_tbl[["pfs_years"]][1] <- 1

  expr_tbl[[2]] <- "some char"
  expect_error(data$qc_preprocess(expr_tbl = expr_tbl), regexp = "numeric")
  expr_tbl[[2]] <- 1

  expr_tbl[[2]][1] <- NA
  expect_error(data$qc_preprocess(expr_tbl = expr_tbl), regexp = "missing values")
  expr_tbl[[2]][1] <- 1

  data$pheno_tbl <- data$pheno_tbl[, c(2:ncol(data$pheno_tbl), 1)]
  expect_error(data$qc_preprocess(expr_tbl = expr_tbl), regexp = "First column")
})
