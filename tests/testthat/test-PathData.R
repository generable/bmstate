test_that("PathData methods work", {
  N <- 100
  pd <- simulate_example_data(N)
  df <- pd$as_data_frame()
  expect_equal(ncol(df), 20)
  df2 <- pd$as_data_frame(truncate = TRUE)
  expect_equal(nrow(df2) + N, nrow(df))

  p1a <- pd$plot_paths()
  p1b <- pd$plot_paths(truncate = TRUE)
  expect_s3_class(p1a, "ggplot")
  expect_s3_class(p1b, "ggplot")
  p2 <- pd$plot_graph()
  p3 <- pd$plot_graph(include_censor = TRUE)

  expect_s3_class(p2, "qgraph")
  expect_s3_class(p3, "qgraph")
  freq <- pd$trans_matrix()
  expect_true(inherits(freq, "matrix"))
})
