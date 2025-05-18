test_that("DosingData init and methods work", {
  N <- 3
  dose_ss <- c(30, 14, 40)
  times <- list(c(4, 5), c(6, 7), c(9, 10))
  doses <- list(c(30, 30), c(0, 30), c(0, 0))
  sid <- c("ASD", "Asd", "EE")
  dd1 <- DosingData$new(sid, doses, times)
  dd2 <- DosingData$new(sid, doses, times, dose_ss, 1)

  theta <- rep(exp(-2), 3)
  theta <- matrix(rep(theta, N), N, 3, byrow = TRUE)
  t <- seq(4, 12, by = 0.1)
  ts <- list()
  for (j in 1:N) {
    ts[[j]] <- t
  }
  a <- dd2$simulate_pk(ts, theta)
  b <- dd2$plot_pk(a)
})
