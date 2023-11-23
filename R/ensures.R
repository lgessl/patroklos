#' @title Ensure patients in expression and pheno data match
#' @description Ensure that the patients in the expression and pheno data match 
#' via subsetting and sorting depenendent on pipeline stage.
#' @param stage string. The pipeline stage we are in. One of `"preprocessing"`,
#' `"after_generate_xy"`.
#' @param expr tibble or numeric matrix.
#' * If `stage == "preprocessing"`: `expr` is a tibble with genes as rows and
#' patients as columns. The column `gene_id_col` (see below) holds the gene
#' dentifiers, all other columns hold the patient identifiers.
#' * If `stage == "after_generate_xy"`: `expr` is a numeric matrix with patients
#' as rows and genes as columns. The row names hold the patient identifiers, the
#' column names hold the gene identifiers.
#' @param pheno tibble or numeric matrix. The pheno data, with patients as rows 
#' and variables as columns. The column `patient_id_col` (see below) holds the
#' patient identifiers.
#' * If `stage == "preprocessing"`: `pheno` is a tibble.
#' * If `stage == "after_generate_xy"`: `pheno` is a numeric matrix as returned
#' by `generate_response()`.
#' @param patient_id_col string. The name of the column in `pheno` that holds
#' the patient identifiers. Default is `NULL`.
#' @param gene_id_col string. The name of the column in `expr` that holds the
#' gene identifiers. `patient_id_col` and `gene_id_col` are only used if
#' `stage == "preprocessing"`. Default is `NULL`.
#' @return A list with two elements, `expr` and `pheno`. Matching and sorted
#' expression and pheno data, i.e. the vectors of patient identifiers are 
#' identical.
#' @export
ensure_patients_match <- function(
    stage,
    expr,
    pheno,
    patient_id_col = NULL,
    gene_id_col = NULL
){
    intersect_ids <- NULL
    cat(nrow(pheno), "samples before matching.\n")
    # case 1: during preprocessing
    if(stage == "preprocessing"){
        patient_ids_expr <- colnames(expr)[-1]
        patient_ids_pheno <- pheno[[patient_id_col]]
        intersect_ids <- intersect(patient_ids_expr, patient_ids_pheno)
        intersect_ids <- sort(intersect_ids)
        # expr: subset and sort
        expr <- expr[, c(gene_id_col, intersect_ids)]
        # pheno: subset and sort
        pheno <- pheno[pheno[[patient_id_col]] %in% intersect_ids, ]
        pheno <- pheno[order(pheno[[patient_id_col]]), ]
    }

    # case 2: after generate_predictor and generate_response
    if("matrix" %in% class(expr) && "matrix" %in% class(pheno)){
        patient_ids_expr <- rownames(expr)
        patient_ids_pheno <- rownames(pheno)
        intersect_ids <- intersect(patient_ids_expr, patient_ids_pheno)
        intersect_ids <- sort(intersect_ids)
        expr <- expr[intersect_ids, ]
        pheno <- pheno[intersect_ids, ]
    }

    if(length(c(patient_ids_expr, patient_ids_pheno)) < 2){
        stop("stage must be one of 'preprocessing' or 'after_generate_xy'.")
    }

    cat(nrow(pheno), "samples after matching.\n")
    res <- list(
        "expr" = expr,
        "pheno" = pheno
    )
    return(res)
}


ensure_available <- function(
    x,
    y
){
    check_consistent_patient_ids(
        stage = "after_generate_xy",
        expr = x,
        pheno = y
    )
    x_fully_available <- apply(x, 1, function(row) any(is.na(row)))
    y_fully_available <- apply(y, 1, function(row) any(is.na(row)))
    x_and_y_fully_available <- x_fully_available & y_fully_available
    x <- x[x_and_y_fully_available, ]
    y <- y[x_and_y_fully_available, ]
    return(list("x" = x, "y" = y))
}