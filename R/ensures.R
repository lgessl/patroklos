#' @title Ensure patients in expression and pheno tibbles match during preprocessing
#' @description Ensure that the patients in the expression and pheno tibbles match 
#' via subsetting and sorting, during preprocessing.
#' @param expr_tbl tibble. With genes as rows and patients as columns. The column 
#' `gene_id_col` (see below) holds the gene identifiers, all other colum names 
#' are the patient identifiers.
#' @param pheno_tbl tibble. The pheno data, with patients as rows and variables as 
#' columns. The column `patient_id_col` (see below) holds the patient identifiers.
#' @param patient_id_col string. The name of the column in `pheno` that holds
#' the patient identifiers.
#' @param gene_id_col string. The name of the column in `expr` that holds the
#' gene identifiers.
#' @param quiet logical. Whether to print messages. Default is `FALSE`.
#' @return A list with two tibbles named `expr` and `pheno`. Matching and sorted
#' expression and pheno data, i.e. the vectors of patient identifiers are 
#' identical.
#' @export
ensure_patients_match <- function(
    expr_tbl,
    pheno_tbl,
    patient_id_col = "patient_id",
    gene_id_col = "gene_id",
    quiet = FALSE
){
    if(!is.data.frame(expr_tbl)){
        stop("expr_tbl must inherit from `data.frame`.")
    }
    if(!is.data.frame(pheno_tbl)){
        stop("pheno_tbl must inherit from `data.frame`.")
    }
    check_tbl_columns_exist(pheno_tbl, "pheno_tbl", patient_id_col)
    check_tbl_columns_exist(expr_tbl, "expr_tbl", gene_id_col)

    if(colnames(expr_tbl)[1] != gene_id_col){
        if(!quiet){
        message("Moving ", gene_id_col, " column to first column")
        }
        expr_tbl <- expr_tbl |> 
            dplyr::relocate(dplyr::all_of(gene_id_col))
    }

    if(!quiet){
        message(nrow(pheno_tbl), " samples in pheno before matching.\n",
            ncol(expr_tbl) - 1, " samples in expr before matching.")
    }
    intersect_ids <- NULL
    patient_ids_expr <- colnames(expr_tbl)[-1]
    patient_ids_pheno <- pheno_tbl[[patient_id_col]]
    intersect_ids <- intersect(patient_ids_expr, patient_ids_pheno)
    intersect_ids <- sort(intersect_ids)
    # expr: subset and sort
    expr_tbl <- expr_tbl[, c(gene_id_col, intersect_ids)]
    # pheno: subset and sort
    pheno_tbl <- pheno_tbl[pheno_tbl[[patient_id_col]] %in% intersect_ids, ]
    pheno_tbl <- pheno_tbl[order(pheno_tbl[[patient_id_col]]), ]

    if(!quiet){
        message(nrow(pheno_tbl), " samples after matching.")
    }
    res <- list(
        "expr_tbl" = expr_tbl,
        "pheno_tbl" = pheno_tbl
    )
    return(res)
}