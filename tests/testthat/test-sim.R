test_that("data simulation works", {
  a <- simulate_example_data(sys_idx = 1)
  b <- simulate_example_data(sys_idx = 2)
  c <- simulate_example_data(sys_idx = 3)
  expect_true(inherits(a$pd, "PathData"))
  expect_true(inherits(b$pd, "PathData"))
  expect_true(inherits(c$pd, "PathData"))
})
