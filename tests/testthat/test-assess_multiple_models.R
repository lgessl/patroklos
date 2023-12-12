test_that("assess_multiple_models() works", {

  set.seed(4352)

  n_samples <- 25
  n_genes <- 5
  n_na_in_pheno <- 5
  n_fold <- 2

  base_dir <- withr::local_tempdir("data")
  data_dir <- file.path(base_dir, "data/mock")
  dir.create(data_dir, recursive = TRUE)
  generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno,
    to_csv = data_dir
  )
  data_spec <- DataSpec(name = "Mock et al. (2023)", directory = data_dir)

  model_dir <- file.path(base_dir, "models")
  dir.create(model_dir)
  model_spec_1 <- ModelSpec(
    name = "cox-zerosum",
    fitter = zeroSum::zeroSum,
    optional_fitter_args = list(family = "cox", alpha = 1, nFold = n_fold, zeroSum = FALSE),
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
    optional_fitter_args = list(family = "binomial", alpha = 1, nFold = n_fold, zeroSum = FALSE),
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

  res_dir <- file.path(base_dir, "results")
  perf_plot_spec <- PerfPlotSpec(
    fname = file.path(res_dir, "perf_plot.pdf"),
    x_metric = "rpp",
    y_metric = "prec",
    show_plots = FALSE
  )

  expect_no_error(
    assess_multiple_models(
        data_spec_list = list(data_spec),
        model_spec_list = model_spec_list,
        perf_plot_spec = perf_plot_spec,
        model_tree_mirror = c("models", "results")
    )
  )
  perf_plot_spec$fname <- file.path(model_dir, "all.pdf")
  expect_no_error(
    assess_train_and_test(
        model_spec_list = list(model_spec_1),
        data_spec_train = data_spec,
        data_spec_test = data_spec,
        perf_plot_spec_train = perf_plot_spec,
        model_tree_mirror = c("models", "results"),
        single_plots = FALSE,
        comparison_plot = FALSE,
        quiet = TRUE
    )
  )
})