test_that("PathData methods work", {
  N <- 100
  pd <- simulate_example_data(N)

  # As data frame
  df1 <- pd$as_data_frame()
  df2 <- pd$as_data_frame(covariates = pd$covariate_names())
  expect_true(ncol(df2) > ncol(df1))
  expect_equal(nrow(df1), nrow(df2))
  df3 <- pd$as_data_frame(truncate = TRUE)
  expect_true(nrow(df3) <= nrow(df1))

  # Plot paths
  p1a <- pd$plot_paths()
  p1b <- pd$plot_paths(truncate = TRUE)

  # As transitions
  dt1 <- pd$as_transitions()
  dt2 <- pd$as_transitions("age")
  expect_true(all(c("from", "to") %in% colnames(dt1)))
  expect_true(ncol(dt1) + 1 == ncol(dt2))

  # As msdata
  dtt <- pd$as_transitions_alt(truncate = TRUE)
  ms <- pd$as_msdata()

  expect_s3_class(p1a, "ggplot")
  expect_s3_class(p1b, "ggplot")
  p2 <- pd$plot_graph()
  p3 <- pd$plot_graph(include_censor = TRUE)

  expect_s3_class(p2, "qgraph")
  expect_s3_class(p3, "qgraph")
  freq <- pd$trans_matrix()
  expect_true(inherits(freq, "matrix"))
})
