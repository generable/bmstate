test_that("MultistateModel init and methods work", {
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
  expect_error(mod$get_knots())

  t_ev <- stats::runif(100)
  mod$set_knots(1, t_ev)
  expect_equal(length(mod$get_knots()), 5)

  w <- rnorm(6)
  w0 <- -4
  ttt <- seq(0, 1, by = 0.01)
  bh1 <- mod$log_baseline_hazard(ttt, w0)
  bh2 <- mod$log_baseline_hazard(ttt, w0, w)
  expect_true(all(bh1 == w0))
  expect_true(length(bh2) == length(ttt))
})


test_that("creating MultistateModel from PathData works", {
  sim <- simulate_example_data(100, sys_idx = 2)
  pd <- sim$pd
  # mod <- create_msm(a, nam, covs, compile = F)
  # expect_true(inherits(mod, "MultistateModel"))
  # expect_output(print(mod))
})
