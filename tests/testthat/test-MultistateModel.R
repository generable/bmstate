test_that("MultistateModel init and methods work", {
  a <- matrix(0, 3, 3)
  a[1, 2] <- 1
  a[2, 3] <- 1
  a[2, 2] <- 1
  states <- c("A", "B", "C")
  tm <- TransitionMatrix$new(a, states)
  covs <- c("CrCL")
  pk_covs <- list(ka = "age", CL = "sex", V2 = "weight")
  mod <- create_msm(tm, covs, pk_covs)
  expect_true(inherits(mod, "MultistateModel"))
  expect_output(print(mod))

  t_ev <- stats::runif(100)
  mod$set_knots(1, t_ev, 9)
  expect_equal(length(mod$system$get_knots()), 9)
  K <- mod$system$num_weights()
  w <- rnorm(K)
  w0 <- -4
  ttt <- seq(0, 1, by = 0.01)
  bh1 <- mod$system$log_baseline_hazard(ttt, w0)
  bh2 <- mod$system$log_baseline_hazard(ttt, w0, w)
  expect_true(all(bh1 == w0))
  expect_true(length(bh2) == length(ttt))
})

test_that("path generation via MultistateSystem works", {
  tm <- transmat_full(self_loops = TRUE)
  mod <- create_msm(tm, "age", NULL)
  tmax <- 3 * 365.25
  mod$set_knots(tmax, seq(0, tmax - 1, length.out = 1000), 4)
  N <- 100
  S <- tm$num_trans()
  L <- mod$system$num_weights()
  w <- array(0.1 * rnorm(N * S * L), dim = c(N, S, L))
  log_w0 <- matrix(-4, N, S)
  log_m <- matrix(0, N, S)
  paths <- mod$system$simulate(w, log_w0, log_m)
  expect_equal(ncol(paths), 6)
})

test_that("data simulation via MultistateModel works", {
  tm <- transmat_full(self_loops = FALSE)
  mod <- create_msm(tm, "age", NULL)
  tmax <- 3 * 365.25
  mod$set_knots(tmax, seq(0, tmax - 1, length.out = 1000), 3)
  N_sub <- 84
  simdat <- mod$simulate_data(N_sub = N_sub)
  expect_true(inherits(simdat, "JointData"))
  expect_equal(simdat$paths$n_paths(), N_sub)
})

test_that("MultistateModel with PK submodel works", {
  tm <- transmat_full(state_names = c("Rand", "Bleed", "Death"))
  covs <- c("sex", "age")
  pk_covs <- list(
    ka = "age",
    CL = "CrCL",
    V2 = c("weight", "sex")
  )
  mod <- create_msm(tm, covs, pk_covs)
  expect_equal(mod$covs(), c("sex", "age", "xpsr"))
  expect_equal(mod$data_covs(), c("sex", "age", "CrCL", "weight"))
  tmax <- 3 * 365.25
  mod$set_knots(tmax, seq(0, tmax - 1, length.out = 1000), 3)
  N <- 99
  ss <- mod$simulate_data(N)
  expect_equal(ss$paths$n_paths(), N)
})
