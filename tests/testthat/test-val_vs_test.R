test_that("val_vs_test() works", { 

    set.seed(123)

    n_samples <- 50
    n_genes <- 2
    
    dir <- withr::local_tempdir()
    data <- generate_mock_data(n_samples = n_samples, n_genes = n_genes)

    model1 <- Model$new(
        name = "regensburg",
        fitter = ptk_zerosum,
        directory = file.path(dir, "model1"),
        time_cutoffs = 2,
        val_error_fun = neg_roc_auc,
        hyperparams = list(family = "gaussian", nFold = 2, lambdaSteps = 2, 
            zeroSum = FALSE)
    )
    model2 <- model1$clone()
    model2$name <- "berlin"
    model2$directory <- file.path(dir, "model2")
    model2$hyperparams[["family"]] <- "binomial"
    model2$time_cutoffs <- 3
    model3 <- model1$clone()
    model3$name <- "frankfurt"
    model3$directory <- file.path(dir, "model3")

    model4 <- Model$new(
        name = "cnn projection",
        directory = file.path(dir, "cnn"),
        fitter = projection_on_feature,
        time_cutoffs = 2,
        val_error_fun = neg_prec_with_prev_greater(0.17),
        hyperparams = list(feature = "ipi"),
        include_from_continuous_pheno = "ipi"
    )

    models <- list(model1, model2, model3, model4)
    data$cohort <- "train"
    training_camp(models, data, quiet = TRUE)

    data$cohort <- "test"
    plt <- val_vs_test(
        models, 
        data, 
        error_fun = neg_roc_auc, 
        regex1 = c("reg", "frank", "."), 
        regex2 = c("ber"), 
        file = file.path(dir, "val_vs_test.pdf"), 
        name1 = c("bavaria", "hesse", "other state"), 
        name2 = c("capital"), 
        quiet = TRUE
    )
    # print(plt)
    expect_s3_class(plt, "gg")
    expect_true(file.exists(file.path(dir, "val_vs_test.pdf")))
    tbl <- val_vs_test(
        models, 
        data, 
        error_fun = neg_prec_with_prev_greater(0.17),
        regex1 = c("reg", "frank"), 
        regex2 = "ber", 
        file = NULL, 
        name1 = c("bavaria", "hesse"), 
        name2 = "capital", 
        correlation_label = FALSE,
        return_type = "tibble", 
        quiet = TRUE
    )
    expect_equal(nrow(tbl), length(models))
    expect_s3_class(tbl, "tbl_df")
    plt <- val_vs_test(
        models,
        data,
        neg_roc_auc,
        regex1 = c("reg", "frank", "."),
        regex2 = NULL,
        file = file.path(dir, "plot.png"),
        name1 = c("bavaria", "hesse", "other state"),
        name2 = NULL,
        legendtitle1 = "federal state",
        legendtitle2 = NULL,
        return_type = "ggplot",
        quiet = TRUE
    )
    print(plt)
    plt <- val_vs_test(
        models,
        data,
        neg_roc_auc,
        regex1 = NULL,
        regex2 = NULL,
        return_type = "ggplot",
        quiet = TRUE
    )
})