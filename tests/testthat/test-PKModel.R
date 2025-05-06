test_that("PKModel init and methods work", {
  covs <- list(ka = "age", CL = "age", V2 = c("weight", "age"))
  a <- PKModel$new(covs)
  expect_message(print(a))

  theta <- list(ka = 0.1, CL = 0.3, V2 = 0.5)
  t <- seq(0, 100, by = 1)
  conc <- a$simulate_ss(t, theta, 100, 24)
  expect_equal(length(conc), length(t))
})
