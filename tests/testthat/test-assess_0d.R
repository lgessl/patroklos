test_that("AssScalar_assess() works", {

    set.seed(342)
    n_samples <- 50
    n_genes <- 5
    n_fold <- 1
    lambda <- 1
    split_index <- 1:2
    
    model_dir <- withr::local_tempdir()
    data <- generate_mock_data(
        n_samples = n_samples,
        n_genes = n_genes,
        n_na_in_pheno = 1,
        split_index = split_index
    )
    expr_mat <- data[["expr_mat"]]
    pheno_tbl <- data[["pheno_tbl"]]

    data <- Data(
        name = "mock",
        directory = "mock",
        train_prop = .66
    )
    model <- Model$new(
        name = "logistic",
        directory = model_dir,
        fitter = zeroSum::zeroSum,
        split_index = split_index,
        time_cutoffs = 2.,
        optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, 
            lambda = lambda, zeroSum = FALSE),
        response_type = "binary"
    )
    ass_scalar <- AssScalar(
        metric = "get_auc",
        pivot_time_cutoff = 2
    )

    prepare_and_fit(expr_mat, pheno_tbl, data, model, quiet = TRUE)

    auc <- AssScalar_assess(
        expr_mat = expr_mat,
        pheno_tbl = pheno_tbl,
        data = data,
        model = model,
        ass_scalar = ass_scalar
    )
    expect_true(is.numeric(auc))
    expect_equal(length(auc), length(split_index))
})