test_that("DosingData init and methods work", {
  N <- 3
  dose_ss <- c(30, 30, 20)
  times <- list(c(4, 5.9), c(6.2, 8.1), c(8.5, 10))
  doses <- list(c(30, 30), c(0, 30), c(0, 0))
  sid <- c("A", "B", "C")
  dd1 <- DosingData$new(sid, doses, times, dose_ss, 2)

  theta <- rep(exp(-2), 3)
  theta <- matrix(rep(theta, N), N, 3, byrow = TRUE)
  t <- seq(0, 11, by = 0.01)
  ts <- list()
  for (j in 1:N) {
    ts[[j]] <- t
  }
  a <- dd1$simulate_pk(ts, theta)
  b <- dd1$plot(a)
  expect_s3_class(b, "ggplot")
})

test_that("DosingData works with zero steady-state doses", {
  N <- 2
  dose_ss <- c(0, 0)
  tm <- seq(0, 120, by = 24)
  times <- list(tm, tm)
  d1 <- rep(30, length(tm))
  d2 <- d1
  d2[5] <- 0
  doses <- list(d1, d2)
  sid <- c("A", "B")
  dd1 <- DosingData$new(sid, doses, times, tau = 24)

  theta <- rep(exp(-2), 3)
  theta <- matrix(rep(theta, N), N, 3, byrow = TRUE)
  t <- seq(0, 150, by = 0.1)
  ts <- list()
  for (j in 1:N) {
    ts[[j]] <- t
  }
  a <- dd1$simulate_pk(ts, theta)
  b <- dd1$plot(a)
  expect_s3_class(b, "ggplot")
})
