test_that("ensure_patients_match function works correctly", {

  expr_tbl <- tibble::tibble(
    gene_id = c("gene1", "gene2"), 
    patient1 = c(1, 2), 
    patient2 = c(3, 4)
  )
  pheno_tbl <- tibble::tibble(
    patient_id = c("patient1", "patient2"), var = c("pheno1", "pheno2"))
  patient_id_col <- "patient_id"
  gene_id_col <- "gene_id"

  expect_message(
    result <- ensure_patients_match(expr_tbl, pheno_tbl, patient_id_col, gene_id_col)
  )
  expect_equal(result[["expr_tbl"]], expr_tbl)
  expect_equal(result[["pheno_tbl"]], pheno_tbl)

  # less benign input
  # need to move gene_id column to first column
  short_expr_tbl <- expr_tbl |> dplyr::relocate(patient2)
  expect_message(
    result <- ensure_patients_match(short_expr_tbl, pheno_tbl, patient_id_col, gene_id_col)
  )
  expect_equal(result[["expr_tbl"]][[1]], expr_tbl[[gene_id_col]])
  # need to subset pheno_tbl
  expr_tbl <- expr_tbl |> dplyr::select(!patient2)
  expect_message(
    result <- ensure_patients_match(expr_tbl, pheno_tbl, patient_id_col, gene_id_col)
  )
  expect_equal(dim(result[["pheno_tbl"]]), c(1L, 2L))
  # need to subset expr_tbl
  short_pheno_tbl <- pheno_tbl |> dplyr::filter(patient_id == "patient1")
  expect_message(
    result <- ensure_patients_match(expr_tbl, short_pheno_tbl, patient_id_col, gene_id_col)
  )
  expect_equal(dim(result[["expr_tbl"]]), c(2L, 2L))
  # wrong input
  expect_error(ensure_patients_match(expr_tbl, pheno_tbl, "wrong_patient_id_col", gene_id_col))
})