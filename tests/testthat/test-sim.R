test_that("data simulation works", {
  possible_covs <- c(
    "sex", "age", "first_dose_amount", "pk_pre_dose", "pk_post_dose",
    "country", "country_num", "dose_adjustment", "dose_arm"
  )
  df_beta_true0 <- data.frame(
    cov_name = rep(c("age", "sex", "first_dose_amount"), each = 4),
    trans_idx = rep(c(1, 2, 3, 4), times = 3),
    beta_true = rep(0, 3 * 4)
  )
  df_beta_true1 <- df_beta_true0
  h0_true <- rep(1e-4, 2)
  df_beta_true1$beta_true[which(df_beta_true1$cov_name == "age")] <- 0.5
  df_beta_true1$beta_true[which(df_beta_true1$cov_name == "sex")] <- 0.5
  a <- simulate_multitransition_data(10, possible_covs, h0_true, df_beta_true1)
  expect_true(inherits(a$pd, "PathData"))
})
