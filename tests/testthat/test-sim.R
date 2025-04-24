test_that("data simulation works", {
  a <- simulate_example_data()
  expect_true(inherits(a$pd, "PathData"))
})
