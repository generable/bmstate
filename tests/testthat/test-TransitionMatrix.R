test_that("TransitionMatrix init and methods work", {
  a <- matrix(0, 3, 3)
  a[1, 2] <- 1
  a[2, 3] <- 1
  a[2, 2] <- 1
  nam <- c("A", "B", "C")
  s <- TransitionMatrix$new(a, nam)
  expect_output(s$print())
  expect_equal(s$source_states(), "A")
  expect_equal(s$absorbing_states(), "C")
  p1 <- s$plot()
  expect_s3_class(p1, "qgraph")
  expect_equal(nrow(s$states_df()), 3)
  expect_equal(nrow(s$trans_df()), 3)
})

test_that("full_transmat() convenience works", {
  a <- full_transmat()
  b <- full_transmat(LETTERS[1:3], self_loops = F)
  expect_equal(nrow(a$trans_df()), 9)
  expect_equal(nrow(b$trans_df()), 3)
})
