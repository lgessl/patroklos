test_that("val_vs_test() works", { 

    set.seed(123)

    n_samples <- 50
    n_genes <- 2
    
    dir <- withr::local_tempdir()
    data <- generate_mock_data(n_samples = n_samples, n_genes = n_genes)

    model1 <- Model$new(
        name = "model1",
        fitter = ptk_zerosum,
        directory = file.path(dir, "model1"),
        time_cutoffs = 2,
        val_error_fun = neg_roc_auc,
        hyperparams = list(family = "gaussian", nFold = 2, lambdaSteps = 2, 
            zeroSum = FALSE)
    )
    model2 <- model1$clone()
    model2$name <- "model2"
    model2$directory <- file.path(dir, "model2")
    model2$hyperparams[["family"]] <- "binomial"
    model2$time_cutoffs <- 3
    model3 <- model1$clone()
    model3$name <- "model3"
    model3$directory <- file.path(dir, "model3")
    models <- list(model1, model2, model3)
    data$cohort <- "train"
    training_camp(models, data, quiet = TRUE)

    data$cohort <- "test"
    expect_no_error(
        plt <- val_vs_test(models, data, error_fun = neg_roc_auc, 
            spotlight_regex = "1|2", file = file.path(dir, "val_vs_test.jpeg"),
            spotlight_name = "leq 2", quiet = TRUE)
    )
    expect_no_error(
        plt <- val_vs_test(models, data, error_fun = neg_roc_auc, 
            spotlight_regex = c("1", "2"), file = file.path(dir, "val_vs_test.jpeg"),
            spotlight_name = c("m1", "m2"), quiet = TRUE)
    )
    print(plt)
})