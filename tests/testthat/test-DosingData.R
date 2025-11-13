test_that("PSSDosingData init and methods work", {
  dose <- c(30, 30, 20)
  sid <- c("A", "B", "C")
  df <- data.frame(subject_id = sid, dose = dose, num_doses = 10, num_ss_doses = 6)
  dd <- simulate_dosing(df, tau = 5, t_jitter = 0.7, p_miss = 0.2)

  theta <- rep(exp(-2), 3)
  N <- length(dose)
  theta <- matrix(rep(theta, N), N, 3, byrow = TRUE)
  t <- seq(0, 60, by = 0.05)
  ts <- list()
  for (j in 1:N) {
    ts[[j]] <- t
  }
  a <- dd$simulate_pk(ts, theta)
  b <- dd$plot(a)
  expect_true(is_ggplot(b))
  r <- dd$filter()
  expect_equal(r$num_subjects(), N)
})

test_that("PSSDosingData works with zero steady-state doses", {
  N <- 1
  dose_ss <- c(0)
  tau <- 12
  tm <- seq(0, 5 * tau, by = tau) + 3
  times <- list(tm)
  d1 <- rep(30, length(tm))
  d1[4] <- 0
  doses <- list(d1)
  sid <- c("A")
  dd1 <- PSSDosingData$new(sid, doses, times, tau = tau)

  theta <- rep(exp(-2), 3)
  theta <- matrix(rep(theta, N), N, 3, byrow = TRUE)
  t <- seq(0, 6 * tau, by = 0.1)
  ts <- list()
  for (j in 1:N) {
    ts[[j]] <- t
  }
  a <- dd1$simulate_pk(ts, theta)
  b <- dd1$plot(a)
  expect_true(is_ggplot(b))
})
