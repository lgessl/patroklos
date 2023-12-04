test_that("read function works correctly", {

  data_spec <- DataSpec(
    name = "Schmitz et al. (2018)",
    directory = "data/schmitz",
    expr_fname = "expr1.csv",
    pheno_fname = "pheno1.csv",
    patient_id_col = "patient_id",
    gene_id_col = "gene_id"
  )

  # Call the function
  result <- read(data_spec)
  expr_mat <- result[["expr_mat"]]
  pheno_tbl <- result[["pheno_tbl"]]

  expect_type(expr_mat, "double")
  expect_true(inherits(expr_mat, "matrix"))
  expect_s3_class(pheno_tbl, "tbl_df")
  expect_equal(colnames(pheno_tbl)[1], data_spec$patient_id_col)

  # errors
  # directory does not exist
  data_spec$directory <- "nonexistent.csv"
  expect_error(read(data_spec))
  data_spec$patient_id_col <- "missing_patient_id_col"
  expect_error(read(data_spec))
  data_spec$gene_id_col <- "missing_gene_id"
  expect_error(read(data_spec))
})