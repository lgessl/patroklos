test_that("discretize_ipi_features() works", {

  set.seed(342)

  n_samples <- 5
  pheno_tbl <- tibble::tibble(
    "age" = rnorm(n_samples, 60, 10),
    "ann_arbor_stage" = sample(1:4, n_samples, replace = TRUE),
    "ldh_ratio" = runif(n_samples, 0, 2),
    "ecog_performance_status" = sample(0:4, n_samples, replace = TRUE),
    "n_extranodal_sites" = sample(0:6, n_samples, replace = TRUE)
  )
  col_names <- colnames(pheno_tbl)
  pheno_tbl[["mock_col"]] <- rnorm(n_samples, 0, 1)
  data <- DataSpec(
    name = "test", 
    directory = "mock", 
    train_prop = 0.8
  )

  new_tbl <- discretize_tbl_cols(
    tbl = pheno_tbl,
    col_names = col_names,
    cutoffs = c(60, 2, 1, 2, 3),
    gl = NULL
  )
  expect_equal(ncol(new_tbl), 5+1+5)
  new_col_bool <- stringr::str_detect(colnames(new_tbl), ">")
  expect_equal(sum(new_col_bool), length(col_names))
  expect_true(all(as.matrix(new_tbl[, new_col_bool]) %in% c(0., 1.)))

  new_tbl <- discretize_tbl_cols(
    tbl = pheno_tbl,
    col_names = col_names,
    cutoffs = c(60, 2, 0, 2, -2),
    gl = c(">", ">", "<", "<", "<")
  )
  expect_equal(ncol(new_tbl), 5+1+5)
  g_bool <- stringr::str_detect(colnames(new_tbl), ">")
  expect_equal(sum(g_bool), 2)
  l_bool <- stringr::str_detect(colnames(new_tbl), "<")
  expect_equal(sum(l_bool), 3)

  expect_error(
    discretize_tbl_cols(
      tbl = pheno_tbl,
      col_names = col_names,
      cutoffs = c(60, 2, 1, 2, 3),
      gl = c(">", ">", "<", "<", "<", ">")
    )
  )
  expect_error(
    discretize_tbl_cols(
      tbl = pheno_tbl,
      col_names = letters[1:5],
      cutoffs = c(60, 2, 1, 2, 3)
    )
  )
})

test_that("prec_from_scores() works", {

  set.seed(345)

  n_samples <- 50
  pheno_tbl <- tibble::tibble(
    "patient" = paste0("patient_", 1:n_samples),
    "IPI" = sample(0:5, n_samples, replace = TRUE),
    "time_to_event" = runif(n_samples, 0, 5),
    "event" = sample(0:1, n_samples, replace = TRUE)
  )
  risk_scores <- runif(n_samples, 0, 1)
  names(risk_scores) <- sample(pheno_tbl[["patient"]])
  names(risk_scores)[1] <- "patient_x"
  data <- DataSpec(
    name = "test", 
    directory = "mock", 
    train_prop = 0.8,
    patient_id_col = "patient",
    time_to_event_col = "time_to_event",
    event_col = "event",
    benchmark_col = "IPI"
  )

  tbl1 <- prec_from_scores(
    pheno_tbl = pheno_tbl,
    data = data,
    risk_scores = risk_scores
  )
  tbl2 <- prec_from_scores(
    pheno_tbl = pheno_tbl,
    data = data
  )
  expect_equal(ncol(tbl1), 3)
  expect_equal(colnames(tbl2), c("IPI >=", "rpp", "prec"))
})

test_that("write_data_info() works", {

  set.seed(345)

  n_samples <- 50
  n_genes <- 50
  pheno_tbl <- tibble::tibble(
    "patient" = paste0("patient_", 1:n_samples),
    "IPI" = sample(0:5, n_samples, replace = TRUE),
    "time_to_event" = runif(n_samples, 0, 5),
    "event" = sample(0:1, n_samples, replace = TRUE)
  )
  risk_scores <- runif(n_samples, 0, 1)
  names(risk_scores) <- sample(pheno_tbl[["patient"]])
  names(risk_scores)[1] <- "patient_x"
  expr_tbl <- tibble::tibble(1:n_genes)
  data <- DataSpec(
    name = "test", 
    directory = "mock", 
    train_prop = 0.8,
    patient_id_col = "patient",
    time_to_event_col = "time_to_event",
    event_col = "event",
    benchmark_col = "IPI"
  )

  filename <- withr::local_tempfile()
  info_json <- write_data_info(
    filename = filename,
    pheno_tbl = pheno_tbl,
    expr_tbl = expr_tbl,
    data = data
  )
  # write into json manually
  info_list <- jsonlite::read_json(filename)
  info_list[["data"]][["expression data"]][["technology"]] <- "rnaseq"
  info_list[["data"]][["benchmark"]][["reference"]] <- "cpp"
  jsonlite::write_json(info_list, filename, auto_unbox = TRUE, pretty = TRUE,
    dataframe = "columns")
  info_json <- write_data_info(
    filename = filename,
    pheno_tbl = pheno_tbl,
    expr_tbl = expr_tbl,
    data = data
  )
  info_list <- jsonlite::read_json(filename)

  expect_equal(info_list[["data"]][["expression data"]][["technology"]][[1]], "rnaseq")
  expect_equal(info_list[["data"]][["benchmark"]][["reference"]][[1]], "cpp")
  expect_equal(length(info_list), 2)
  expect_false(is.null(info_list[["data"]]))
  expect_false(is.null(info_list[["data"]][["pheno data"]]))
  expect_false(is.null(info_list[["data"]][["expression data"]]))
  expect_false(is.null(info_list[["data"]][["benchmark"]]))
  expect_false(is.null(info_list[["data"]][["benchmark"]][["performance"]]))
})
