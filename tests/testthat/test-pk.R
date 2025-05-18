test_that("exposed PK functions work", {
  N <- 3
  dose_ss <- c(30, 14, 40)
  times <- list(c(4, 5), c(6, 7), c(9, 10))
  doses <- list(c(30, 30), c(0, 30), c(0, 0))
  theta <- rep(1, 3)
  theta <- matrix(rep(theta, N), N, 3, byrow = TRUE)
  a <- pop_2cpt_partly_ss(NULL, dose_ss, times, doses, theta)
})
