test_that("assess_model() works", {
  
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
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(dir, "rpp.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    pivot_time_cutoff = 2.,
    benchmark = "ipi",
    show_plots = TRUE
  )

  for(model_spec in list(model_spec_1, model_spec_2)){
    prepare_and_fit(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data_spec = data_spec,
      model_spec = model_spec
    )
  }

  tbl <- assess_model(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_1,
    perf_plot_spec = perf_plot_spec
  )$data
  expect_s3_class(tbl, "tbl_df")

  
  perf_plot_spec$benchmark <- "ipi"
  perf_plot_spec$directory <- file.path(dir, "logrank.pdf")
  perf_plot_spec$y_metric <- "logrank"
  perf_plot_spec$scale_y <- "log10"
  tbl <- assess_model(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data_spec = data_spec,
    model_spec = model_spec_2,
    perf_plot_spec = perf_plot_spec
  )$data
  expect_equal(names(tbl), c("rpp", "logrank", "cutoff", "split", "model"))
})
