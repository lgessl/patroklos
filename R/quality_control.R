qc_prepare <- function(
    x,
    y
){
    # check class and type of x and y
    x_y_list <- list("x" = x, "y" = y)
    for(x_y in names(x_y_list)){
        if(!(is.matrix(x_y_list[[x_y]]))){
            stop(x_y, "must be a matrix.")
        }
        if(!is.numeric(x_y_list[[x_y]])){
            stop(x_y, " must be numeric.")
        }
    }

    check_consistent_patient_ids(
        stage = "after_generate_xy",
        expr = x,
        pheno = y
    )
    check_available(
        x = x,
        y = y
    )
}


qc_preprocess <- function(
    expr_tbl,
    pheno_tbl,
    data_spec,
    check_default = FALSE
){
    if(check_default){
        directory <- data_spec$directory
        data_spec <- DataSpec(
            name = "default",
            directory = directory
        )
    }

    # Extract
    directory <- data_spec$directory
    expr_fname <- data_spec$expr_fname
    pheno_fname <- data_spec$pheno_fname
    gene_id_col <- data_spec$gene_id_col
    patient_id_col <- data_spec$patient_id_col
    pfs_col <- data_spec$pfs_col
    progression_col <- data_spec$progression_col
    ipi_col <- data_spec$ipi_col

    # Check if files exist
    expr_fname <- file.path(directory, expr_fname)
    pheno_fname <- file.path(directory, pheno_fname)
    if(!file.exists(expr_fname)){
        stop("Expression file ", expr_fname, " does not exist.")
    }
    if(!file.exists(pheno_fname)){
        stop("Pheno file ", pheno_fname, " does not exist.")
    }

    # Expression tibble
    if(names(expr_tbl)[1] != gene_id_col){
        stop("First column of expression tibble must be ", gene_id_col)
    }
    if(!is.character(expr_tbl[[gene_id_col]])){
        stop("Gene ids must be characters.")
    }
    if(!all(sapply(expr_tbl[, -1], is.numeric))){
        stop("Expression values must be numeric.")
    }
    if(any(is.na(expr_tbl[, -1]))){
        stop("Expression values must not contain missing values.")
    }

    # Pheno tibble
    check_tbl_columns_exist(
        tbl = pheno_tbl,
        tbl_name = "pheno_tbl",
        col_names = c(
            patient_id_col,
            pfs_col,
            progression_col,
            ipi_col
        )
    )
    if(!is.character(pheno_tbl[[patient_id_col]]) || !elements_unique(pheno_tbl[[patient_id_col]])){
        stop("Patient ids must be unique characters.")
    }
    if(!is.numeric(pheno_tbl[[pfs_col]])){
        stop("PFS values must be numeric and not contain missing values.")
    }
    if(any(is.na(pheno_tbl[[pfs_col]]))){
        warning("PFS values contain missing values.")
    }
    if(!is.numeric(pheno_tbl[[progression_col]]) || 
        !all(pheno_tbl[[progression_col]] %in% c(0, 1))){
        stop("Progression values must be numeric and either 1 (progression) or 0 (no progression).")
    }
    if(!is.numeric(pheno_tbl[[ipi_col]]) || !all(pheno_tbl[[ipi_col]] %in% c(0:5, NA))){
        stop("IPI values must be numeric values in 0:5")
    }

    # Both
    check_consistent_patient_ids(
        stage = "preprocessing",
        expr = expr_tbl,
        pheno = pheno_tbl,
        data_spec = data_spec
    )
}