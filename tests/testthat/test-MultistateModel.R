test_that("MultistateModel init works", {
  a <- matrix(0, 3, 3)
  a[1, 2] <- 1
  a[2, 3] <- 1
  a[2, 2] <- 1
  nam <- c("A", "B", "C")
  covs <- list(hazard = "age")
  mod <- create_msm(a, nam, covs, compile = F)
  expect_true(inherits(mod, "MultistateModel"))
  expect_message(print(mod))
})


test_that("creating MultistateModel from PathData works", {
  sim <- simulate_example_data(100, sys_idx = 2)
  pd <- sim$pd
  mod <- create_msm(a, nam, covs, compile = F)
  expect_true(inherits(mod, "MultistateModel"))
  expect_message(print(mod))
})
