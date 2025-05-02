test_that("StateGraph init and methods work", {
  a <- matrix(0, 3, 3)
  a[1, 2] <- 1
  a[2, 3] <- 1
  a[2, 2] <- 1
  nam <- c("A", "B", "C")
  s <- StateGraph$new(a, nam)
})
