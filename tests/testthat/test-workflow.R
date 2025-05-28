# Options
options <- list(
  N_subject = 100,
  iter_warmup = 60,
  iter_sampling = 30,
  chains = 1
)

test_that("entire workflow works", {
  # Setup
  h0_base <- 1e-3
  setup <- example_sim_setup_illnessdeath()
  mod <- setup$model

  # Simulate data
  jd <- mod$simulate_data(options$N_subject, beta_haz = setup$beta_haz)

  # Split
  jd <- split_data(jd)
  expect_true(inherits(jd$test, "JointData"))
  expect_equal(jd$test$paths$n_paths(), options$N_subject * 0.25)

  # CoxPH fit
  cph <- jd$train$paths$fit_coxph(covariates = mod$covs())
  msf <- jd$train$paths$fit_mstate()
  expect_true(inherits(msf, "msfit"))

  # Fit the model
  a <- fit_stan(mod, jd$train,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains,
    refresh = 5,
    adapt_delta = 0.95,
    init = 0.1,
    return_stanfit = TRUE
  )
  fit <- a$fit
  expect_true(inherits(fit, "MultistateModelStanFit"))

  # Computing hazard multipliers
  log_m <- msmsf_log_hazard_multipliers(fit)
  log_m_test <- msmsf_log_hazard_multipliers(fit, data = jd$test)

  # Path generation
  p <- generate_paths(fit, n_rep = 3)
  p_oos <- generate_paths(fit, n_rep = 3, data = jd$test)

  pe <- p_event(p)
  pe_bysex <- p_event(p, by = "sex")
  pe_bydose <- p_event(p, by = "dose_amt")
  expect_equal(nrow(pe), 3)
  expect_equal(nrow(pe_bysex), 5)
  expect_equal(nrow(pe_bydose), 7)

  # Baseline hazard viz
  plot_h0 <- fit$plot_h0() # should be at h0_base level
  expect_s3_class(plot_h0, "ggplot")
})



test_that("entire workflow works (with PK)", {
  # Setup
  h0_base <- 1e-3

  # Simulate data
  sim <- simulate_example_data(N = options$N_subject)
  mod <- sim$model
  jd <- sim$data

  # Split
  jd <- split_data(jd)

  # Fit the model
  fit <- fit_stan(mod, jd$train,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains,
    refresh = 5,
    adapt_delta = 0.95,
    init = 0.1
  )
  expect_true(inherits(fit, "MultistateModelStanFit"))

  # Pk params
  pkpar <- msmsf_pk_params(fit)
  pkpar_oos <- msmsf_pk_params(fit, jd$test)
  expect_true(length(pkpar) == 30)
  expect_true(length(pkpar_oos) == 30)

  # Hazard multipliers
  log_m <- msmsf_log_hazard_multipliers(fit)
  log_m_test <- msmsf_log_hazard_multipliers(fit, jd$test)
  N_sub_test <- length(jd$test$paths$unique_subjects())
  H <- mod$system$num_trans()
  expect_equal(dim(log_m_test[[1]]), c(N_sub_test, H))

  # Path prediction
  p <- generate_paths(fit)
  P <- solve_trans_prob_fit(fit)

  # Paths from 0 to 3 years starting from Randomization

  # Paths from 0 to 3 years starting from Randomization (oos)

  # Covariate effect and baseline hazard plots
  # plt_oth <- plot_other_beta(res$fit, res$stan_dat, res$jd, res$df_beta_true)

  # plt_h0 <- plot_h0(res$fit, sd, res$jd, df_h0)
  # expect_s3_class(plt_h0, "ggplot")

  # Event probabilities
  # names <- c("1. Obs", "2. IS pred", "3. OOS pred")
  # plt_a <- plot_p_event(res$jd, res$paths_3yr, res$paths_3yr_oos, names)

  # no points selected for one or more curves, consider using the extend argument
  # plt_b <- plot_p_event_by_dose(res$jd, res$paths_3yr, res$paths_3yr_oos, names)
  # expect_s3_class(plt_a, "ggplot")
  # expect_s3_class(plt_b, "ggplot")

  # Event summary
  # There was 1 warning in `dplyr::mutate()`.
  # i In argument: `event_surv = furrr::future_map(...)`.
  # Caused by warning in `serializedSize()`:
  #  ! 'package:stats' may not be available when loading
  # target_times <- seq(0, 3 * 365.25, by = 365.25 / 4)
  # ppsurv_subj <- summarize_ppsurv(
  #  res$paths_3yr,
  #  target_times = target_times,
  #  by = c("subject_id"),
  #  truncate_at_terminal_events = FALSE
  # )
  # expect_true(inherits(ppsurv_subj, "tbl_df"))
})
