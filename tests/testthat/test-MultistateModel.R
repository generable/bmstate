test_that("MultistateModel init works", {
  a <- matrix(0, 3, 3)
  a[1, 2] <- 1
  a[2, 3] <- 1
  a[2, 2] <- 1
  nam <- c("A", "B", "C")
  mod <- create_msm(a, nam, compile = F)
  expect_true(inherits(mod, "MultistateModel"))
  expect_output(mod$print())
})
