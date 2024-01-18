test_that("read function works correctly", {

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
  data_spec <- DataSpec(
    name = "mock",
    directory = dir,
    train_prop = .7,
    expr_fname = "expr.csv",
    pheno_fname = "pheno.csv",
    patient_id_col = "patient_id",
    gene_id_col = "gene_id"
  )

  # Call the function
  result <- read(data_spec)
  expr_mat <- result[["expr_mat"]]
  pheno_tbl <- result[["pheno_tbl"]]

  expect_true(is.numeric(expr_mat))
  expect_true(is.matrix(expr_mat))
  expect_s3_class(pheno_tbl, "tbl_df")
  expect_equal(colnames(pheno_tbl)[1], data_spec$patient_id_col)

  # errors
  # directory does not exist
  data_spec$patient_id_col <- "missing_patient_id_col"
  expect_error(read(data_spec), "no column named")
  data_spec$gene_id_col <- "missing_gene_id"
  expect_error(read(data_spec), "no column named")
})
