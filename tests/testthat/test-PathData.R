test_that("PathData methods work", {
  N <- 100
  pd <- simulate_example_data(N)
  expect_output(print(pd))

  # As data frame
  df1 <- pd$as_data_frame()
  df2 <- pd$as_data_frame(covariates = pd$covariate_names())
  expect_true(ncol(df2) > ncol(df1))
  expect_equal(nrow(df1), nrow(df2))
  df3 <- pd$as_data_frame(truncate = TRUE)
  expect_true(nrow(df3) <= nrow(df1))

  # As transitions
  dt1 <- pd$as_transitions()
  dt2 <- pd$as_transitions("age")
  expect_true(all(c("from", "to") %in% colnames(dt1)))
  expect_true(ncol(dt1) + 1 == ncol(dt2))

  # As msdata
  dtt <- pd$as_transitions_alt(truncate = TRUE)
  ms <- pd$as_msdata()
  expect_true(inherits(ms, "msdata"))

  # Plots
  p1a <- pd$plot_paths()
  p1b <- pd$plot_paths(truncate = TRUE)
  expect_s3_class(p1a, "ggplot")
  expect_s3_class(p1b, "ggplot")
  p2 <- pd$plot_graph()
  expect_s3_class(p2, "qgraph")
  prop <- pd$prop_matrix()
  expect_true(inherits(prop, "table"))

  # Cox PH fit
  cp <- pd$coxph(c("age", "sex"))
  expect_true(inherits(cp, "coxph"))
})
