test_that("columns check works", {
  df <- data.frame(hi = c(3, 4))
  expect_true(check_columns(df, "hi"))
  expect_error(check_columns(df, c("asd")))
})
