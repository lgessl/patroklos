test_that("prepare_and_fit", {

  set.seed(4352)

  n_samples <- 25
  n_genes <- 5
  n_na_in_pheno <- 5
  n_fold <- 2

  data_dir <- withr::local_tempdir()
  generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = data_dir
  )
  data_spec <- DataSpec(name = "Mock et al. (2023)", directory = data_dir)

  model_dir <- withr::local_tempdir()
  model_spec_1 <- ModelSpec(
    name = "cox-zerosum",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold),
    response_type = "survival_censored",
    include_from_continuous_pheno = NULL,
    include_from_discrete_pheno = NULL,
    save_dir = file.path(model_dir, "cox"),
    pfs_leq = 2.,
    fit_fname = "model1.rds"
  )
  model_spec_2 <- ModelSpec(
    name = "binomial-zerosum",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold),
    response_type = "binary",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var",
    save_dir = file.path(model_dir, "binomial"),
    pfs_leq = 2.,
    fit_fname = "model2.rds"
  )
  model_spec_list <- list(model_spec_1, model_spec_2)

  data <- read(data_spec)
  prepare_and_fit(
    expr_mat = data[["expr_mat"]],
    pheno_tbl = data[["pheno_tbl"]],
    data_spec = data_spec,
    model_spec_list = model_spec_list
  )

  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(model_dir, "perf_plot.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    show_plots = TRUE
  )

  expect_no_error(
    compare_models(
        data_spec_list = list(data_spec),
        model_spec_list = model_spec_list,
        perf_plot_spec = perf_plot_spec
    )
  )
})