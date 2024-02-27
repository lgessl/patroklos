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
  pheno_tbl <- data$pheno_tbl
  expr_tbl <- data$expr_mat |> t() |> tibble::as_tibble(rownames = "gene_id")
  data <- Data$new(
    name = "mock",
    directory = dir,
    train_prop = .7,
    
  )

  qc_preprocess(
    data = data,
    expr_tbl = expr_tbl,
    pheno_tbl = pheno_tbl
  )
  qc_preprocess(
    data = data,
    expr_tbl = expr_tbl,
    pheno_tbl = pheno_tbl
  )

  data$expr_file <- "nonexistent.csv"
  expect_error(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "Expression file"
  )
  data$expr_file <- "expr.csv"

  pheno_tbl[["progression"]][1] <- 2
  expect_error(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "either 1"
  )
  pheno_tbl[["progression"]][1] <- 0

  pheno_tbl[["patient_id"]][1] <- "sample_2"
  expect_error(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "Patient ids"
  )
  pheno_tbl[["patient_id"]][1] <- "sample_1"

  pheno_tbl[["pfs_years"]][1] <- NA
  expect_warning(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "missing values"
  )
  pheno_tbl[["pfs_years"]][1] <- 1

  expr_tbl[[2]] <- "some char"
  expect_error(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "numeric"
  )
  expr_tbl[[2]] <- 1

  expr_tbl[[2]][1] <- NA
  expect_error(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "missing values"
  )
  expr_tbl[[2]][1] <- 1

  pheno_tbl <- pheno_tbl[, c(2:ncol(pheno_tbl), 1)]
  expect_error(
    qc_preprocess(
      expr_tbl = expr_tbl,
      pheno_tbl = pheno_tbl,
      data = data
    ),
    regexp = "First column"
  )
})
