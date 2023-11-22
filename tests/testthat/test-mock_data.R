test_that("mock_data_from_existing returns correct output", {

    directory <- test_path("data")
    n_samples <- 5L
    n_genes <- 3L
    expr_fname <- "expr.csv"
    pheno_fname <- "pheno.csv"
    patient_id_col <- "patient_id"
    gene_id_col <- "gene_id"
    save <- TRUE
    save_suffix <- "mock"
    cleanup <- TRUE

    # generate mock real data (if not there yet)
    if(!file.exists(file.path(directory, expr_fname)) ||
        !file.exists(file.path(directory, pheno_fname))){    
        generate_mock_data(
            directory = directory,
            n_samples = 2*n_samples,
            n_genes = 2*n_genes,
            expr_fname = expr_fname,
            pheno_fname = pheno_fname,
            patient_id_col = patient_id_col,
            gene_id_col = gene_id_col
        )
    }

    # Call the function
    result <- mock_data_from_existing(
        directory, 
        n_samples, 
        n_genes, 
        expr_fname, 
        pheno_fname, 
        patient_id_col, 
        gene_id_col, 
        save, 
        save_suffix
    )
    expr_tbl <- result[["expr"]]
    pheno_tbl <- result[["pheno"]]

    # tests
    expect_identical(length(result), 3L)
    expect_identical(dim(expr_tbl), c(n_genes, n_samples+1L))
    expect_equal(dim(pheno_tbl)[1], n_samples)
    expect_true(all(
        colnames(expr_tbl)[2:ncol(expr_tbl)] == pheno_tbl[[patient_id_col]]
    ))
    expect_identical(result[["file_names"]], c("expr_mock.csv", "pheno_mock.csv"))

    # clean up
    if(cleanup){
        for(i in 1:2){
            file <- file.path(directory, result[["file_names"]][i])
            file.remove(file)
        }
    }
})
