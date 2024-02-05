test_that("assess_2d() works", {
  
  set.seed(132)
  n <- 20

  dir <- withr::local_tempdir()
  n_fold <- 3
  lambda <- 1
  data <- generate_mock_data(
    n_samples = n,
    n_genes = 5,
    n_na_in_pheno = 1
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]
  pheno_tbl[["split_1"]] <- sample(
    c("train", "test"),
    size = n,
    replace = TRUE
  )
  pheno_tbl[["split_2"]] <- pheno_tbl[["split_1"]]

  data_spec <- DataSpec(
    name = "Mock et al. (2023)",
    directory = "mock",
    train_prop = .66
  )
  model_spec_1 <- ModelSpec(
    name = "cox",
    directory = file.path(dir, "cox"),
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, 
      lambda = lambda, zeroSum = FALSE),
    response_type = "survival_censored"
  )
  model_spec_2 <- ModelSpec(
    name = "logistic",
    directory = file.path(dir, "logistic"),
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2.,
    optional_fitter_args = list(family = "binomial", alpha = 1, 
      nFold = n_fold, lambda = lambda, zeroSum = FALSE),
    response_type = "binary"
  )
  ass_2d_spec <- Ass2dSpec(
    fname = file.path(dir, "rpp.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    pivot_time_cutoff = 2.,
    benchmark = "ipi",
    show_plots = FALSE,
    text = list(ggplot2::aes(x = .5, y = .5, label = "this text"), 
      color = "red", angle = 45)
  )

  for(model_spec in list(model_spec_1, model_spec_2)){
    prepare_and_fit(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
    )
  }

  tbl <- assess_2d(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_1,
    ass_2d_spec = ass_2d_spec
  )$data
  expect_s3_class(tbl, "tbl_df")

  
  ass_2d_spec$benchmark <- "ipi"
  ass_2d_spec$directory <- file.path(dir, "logrank.pdf")
  ass_2d_spec$y_metric <- "logrank"
  ass_2d_spec$scale_y <- "log10"
  tbl <- assess_2d(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_2,
    ass_2d_spec = ass_2d_spec
  )$data
  expect_equal(names(tbl), c("rpp", "logrank", "cutoff", "split", "model"))

  ass_2d_spec$y_metric <- "precision_ci"
  ass_2d_spec$ci_level <- .95
  ass_2d_spec$directory <- file.path(dir, "precision_ci.pdf")
  ass_2d_spec$scale_y <- "identity"
  ass_2d_spec$title <- "Lower precision CI boundary (upper for ipi)"
  ass_2d_spec$show_plots <- FALSE
  tbl <- assess_2d(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_2,
    ass_2d_spec = ass_2d_spec
  )$data
  expect_equal(names(tbl), c("rpp", "precision_ci", "cutoff", "split", "model"))
})
