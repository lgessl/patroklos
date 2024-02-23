test_that("prepare() works",{

  set.seed(4352)

  n_samples <- 30
  n_genes <- 5
  n_na_in_pheno <- 5

  data <- generate_mock_data(
    n_samples = n_samples,
    n_genes = n_genes,
    n_na_in_pheno = n_na_in_pheno
  )
  expr_mat <- data[["expr_mat"]]
  pheno_tbl <- data[["pheno_tbl"]]
  data <- DataSpec(
    name = "Mock et al. (2023)",
    directory = "some_dir",
    train_prop = 0.8,
    cohort = "train"
  )
  model <- ModelSpec(
    name = "zerosum",
    directory = "some_dir",
    fitter = zeroSum::zeroSum,
    split_index = 1,
    time_cutoffs = 2.,
    response_type = "survival_censored",
    include_from_continuous_pheno = "continuous_var",
    include_from_discrete_pheno = "discrete_var"
  )

  expect_no_error(
    result <- prepare(
      expr_mat = expr_mat,
      pheno_tbl = pheno_tbl,
      data = data,
      model = model
    )
  )

  model$response_type <- "binary"
  model$include_from_continuous_pheno <- NULL
  result <- prepare(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data = data,
    model = model
  )

  colnames(pheno_tbl)[1] <- "patient"
  colnames(pheno_tbl)[3] <- "pfs"
  data$patient_id_col <- "patient"
  data$time_to_event_col <- "pfs"
  model$pivot_time_cutoff <- 2.3
  result <- prepare(
    expr_mat = expr_mat,
    pheno_tbl = pheno_tbl,
    data = data,
    model = model
  )
})
