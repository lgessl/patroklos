#' @title Read expression and pheno data from csv files into consistent format
#' @description Read expression data into a matrix and pheno data into a tibble.
#' Both will hold patients as rows.
#' @param data_spec DataSpec S3 object. Specifications on the data. See the the
#' constructor `DataSpec()` for details.
#' @return A list with a numeric matrix, named `expr`, and a tibble named `pheno`. 
#' `expr` holds the expression data, with patient ids as row names and gene ids as
#' column names. I.e., we transpose the expression data. `pheno` holds the pheno data, 
#' with the patient ids in the first column `patient_id_col`.
read <- function(
    data_spec
){
    if(!inherits(data_spec, "DataSpec")){
        stop("data_spec must be a DataSpec object")
    }
    
    # extract values from data_spec
    directory <- data_spec$directory
    expr_file <- data_spec$expr_file
    pheno_file <- data_spec$pheno_file
    patient_id_col <- data_spec$patient_id_col
    gene_id_col <- data_spec$gene_id_col

    # read
    files <- c(expr_file, pheno_file)
    tbls <- list()
    for (i in 1:length(files)){
        full_path <- file.path(directory, files[i])
        tbls[[i]] <- readr::read_csv(full_path, show_col_types = FALSE)
    }
    expr_tbl <- tbls[[1]]
    pheno_tbl <- tbls[[2]]

    # check if identifier columns are there
    if(is.null(expr_tbl[[gene_id_col]])){
        stop("There is no column named ", gene_id_col, " in ", expr_file)
    }
    if(is.null(pheno_tbl[[patient_id_col]])){
        stop("There is no column named ", patient_id_col, " in ", pheno_file)
    }

    # expression to matrix
    gene_names <- tbls[["expr"]][[gene_id_col]]
    if(!elements_unique(gene_names)){
        stop("Column ", gene_id_col, " in ", expr_file, " holds duplicate entries.")
    }
    expr_mat <- expr_tbl |>
        dplyr::select(!dplyr::all_of(gene_id_col)) |>
        as.matrix() |> 
        t()
    colnames(expr_mat) <- gene_names

    # pheno: move patient ids into first column
    patient_ids <- pheno_tbl[[patient_id_col]]
    if(!elements_unique(patient_ids)){
        stop("Column ", patient_id_col, " in ", pheno_file, " holds duplicate entries.")
    }
    pheno_tbl <- pheno_tbl |>
        dplyr::relocate(dplyr::all_of(patient_id_col))

    res <- list(
        "expr_mat" = expr_mat, 
        "pheno_tbl" = pheno_tbl
        )
    return(res)
}