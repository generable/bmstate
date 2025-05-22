test_that("PKModel init and methods work", {
  covs <- list(ka = "age", CL = "age", V2 = c("weight", "age"))
  a <- PKModel$new(covs)
  expect_output(print(a))

  theta <- list(ka = 0.1, CL = 0.3, V2 = 0.5)
  t <- seq(0, 100, by = 1)
  conc <- a$simulate_ss(t, theta, 100, 24)
  expect_equal(length(conc), length(t))
  r <- a$format_params()
  expect_equal(names(r), c("ka", "CL", "V2"))
  expect_equal(a$covs(), c("age", "weight"))
})
