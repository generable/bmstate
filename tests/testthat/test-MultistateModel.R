test_that("MultistateModel init works", {
  a <- matrix(0, 3, 3)
  a[1, 2] <- 1
  a[2, 3] <- 1
  a[2, 2] <- 1
  states <- c("A", "B", "C")
  covs <- c("CrCL")
  pk_covs <- list(ka = "age", CL = "sex", V2 = "weight")
  mod <- create_msm(a, states, covs, pk_covs, compile = F)
  expect_true(inherits(mod, "MultistateModel"))
  expect_output(print(mod))
})


test_that("creating MultistateModel from PathData works", {
  sim <- simulate_example_data(100, sys_idx = 2)
  pd <- sim$pd
  # mod <- create_msm(a, nam, covs, compile = F)
  # expect_true(inherits(mod, "MultistateModel"))
  # expect_output(print(mod))
})
