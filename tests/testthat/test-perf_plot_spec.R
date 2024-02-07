test_that("infer_as2() works", {

  mirror <- c("models", "results")
  ass_spec_2d <- AssSpec2d(
    file = "ass_spec_2d.png",
    x_metric = "x",
    y_metric = "y",
    pivot_time_cutoff = 2.,
    benchmark = "benchmark",
    fellow_csv = TRUE
  )
  ass_spec_2d$model_tree_mirror <- mirror
  data_spec <- DataSpec(
    name = "mock",
    directory = "data/mock",
    train_prop = .7,
    cohort = "train"
  )
  model_spec <- ModelSpec(
    name = "cox",
    directory = "models/cox",
    fitter = zeroSum::zeroSum,
    split_index = 1:2,
    time_cutoffs = 2.
  )

  auto_pps <- infer_as2(
    ass_spec_2d = ass_spec_2d,
    model_spec = model_spec,
    data_spec = data_spec
  )
  expect_equal(
    auto_pps$file,
    "models/cox/x_vs_y.png"
  )
  expect_equal(
    auto_pps$title,
    "mock train, pfs_years < 2"
  )

  data_spec$cohort <- "test"
  auto_pps <- infer_as2(
    ass_spec_2d = ass_spec_2d,
    model_spec = model_spec,
    data_spec = data_spec
  )
  expect_equal(
    auto_pps$file,
    "results/cox/x_vs_y.png"
  )
  expect_equal(
    auto_pps$title,
    "mock test, pfs_years < 2"
  )
})
