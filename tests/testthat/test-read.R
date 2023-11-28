test_that("read function works correctly", {

  directory <- test_path("data/schmitz")
  expr_fname <- "expr.csv"
  pheno_fname <- "pheno.csv"
  patient_id_col <- "patient_id"
  gene_id_col <- "gene_id"
  data_spec <- DataSpec(
    patient_id_col = patient_id_col,
    gene_id_col = gene_id_col
  )

  # Call the function
  result <- read(
    directory, 
    expr_fname, 
    pheno_fname, 
    data_spec
    )
  expr_mat <- result[["expr_mat"]]
  pheno_tbl <- result[["pheno_tbl"]]

  expect_type(expr_mat, "double")
  expect_true(inherits(expr_mat, "matrix"))
  expect_s3_class(pheno_tbl, "tbl_df")
  expect_equal(colnames(pheno_tbl)[1], patient_id_col)

  # errors
  expect_error(read("nonexistent.csv", expr_fname, pheno_fname, patient_id_col, gene_id_col))
  data_spec$patient_id_col <- "missing_patient_id_col"
  expect_error(read(directory, expr_fname, pheno_fname, data_spec))
  data_spec$missing_gene_id_col <- "gene_id"
  expect_error(read(directory, expr_fname, pheno_fname, data_spec))
})