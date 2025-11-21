test_that("JointData methods work", {
  N <- 100
  simdat <- simulate_example_data(N)
  jd <- simdat$data
  expect_true(inherits(jd, "JointData"))
  expect_output(print(jd))
})


test_that("PathData methods work", {
  N <- 100
  simdat <- simulate_example_data(N)
  pd <- simdat$data$paths
  expect_output(print(pd))

  # Covariate subset
  pd2 <- pd$subset_covariates("weight")
  expect_equal(pd2$covs, "weight")

  # As any event
  pd3 <- as_any_event(pd, null_state = "Alive")
  expect_true("Alive" %in% pd3$state_names())

  # As data frame
  df1 <- pd$as_data_frame()
  df2 <- pd$as_data_frame(covariates = pd$covariate_names())
  expect_true(ncol(df2) > ncol(df1))
  expect_equal(nrow(df1), nrow(df2))
  df3 <- pd$as_data_frame(truncate = TRUE)
  expect_true(nrow(df3) <= nrow(df1))

  # As transitions
  dt1 <- pd$as_transitions()
  dt2 <- pd$as_transitions("age")
  expect_true(all(c("from", "to") %in% colnames(dt1)))
  expect_true(ncol(dt1) + 1 == ncol(dt2))

  # As msdata
  dtt <- pd$as_transitions_alt(truncate = TRUE)
  ms <- pd$as_msdata()
  expect_true(inherits(ms, "msdata"))

  # Plots
  p1a <- pd$plot_paths()
  p1b <- pd$plot_paths(truncate = TRUE)
  expect_true(is_ggplot(p1a))
  expect_true(is_ggplot(p1b))
  p2 <- pd$plot_graph()
  prop <- pd$prop_matrix()
  expect_true(inherits(prop, "table"))

  # Frequentist fits
  cp <- pd$fit_coxph(c("age", "sex"))
  expect_true(inherits(cp, "coxph"))
  msf <- pd$fit_mstate()
  leg <- pd$transmat$trans_df()
  plt <- msfit_plot_cumhaz(msf, leg)
  expect_true(is_ggplot(plt))
})
