test_that("PathData methods work", {
  a <- simulate_example_data(100)
  pd <- a$pd
  plt <- pd$plot_paths()
  expect_s3_class(plt, "ggplot")
  freq <- pd$trans_matrix()
  expect_true(inherits(freq, "matrix"))
  expect_null(pd$plot_graph())
})
