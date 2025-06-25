# Options
options <- list(
  N_subject = 100,
  iter_warmup = 100,
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
  expect_true(inherits(fit, "MultistateModelFit"))

  # Computing hazard multipliers
  log_m <- msmsf_log_hazard_multipliers(fit)
  log_m_test <- msmsf_log_hazard_multipliers(fit, data = jd$test)

  # Path generation
  p <- generate_paths(fit, n_rep = 3)
  p_oos <- generate_paths(fit, n_rep = 3, data = jd$test)

  pe <- p_event(p)
  pe_oos <- p_event(p_oos)
  pe_bysex <- p_event(p, by = "sex")
  pe_bydose <- p_event(p, by = "dose_amt")
  expect_equal(nrow(pe), 2)
  expect_equal(nrow(pe_bysex), 4)
  expect_equal(nrow(pe_bydose), 6)

  # Baseline hazard viz
  plot_h0 <- fit$plot_h0() # should be at h0_base level
  expect_s3_class(plot_h0, "ggplot")

  # Scoring
  ev <- "Death"
  er <- event_risk(p, ev)
  a <- as_survival(jd$train$paths, ev)
  a <- a |> dplyr::left_join(er, by = "subject_id")
  ci <- c_index(a)
  expect_gt(ci$concordance, 0)

  # Test that reducing to mean works
  mfit <- fit$mean_fit()
  expect_equal(mfit$num_draws(), 1)
  p_mfit <- generate_paths(mfit, n_rep = 100)
  expect_equal(p_mfit$n_paths(), 100 * jd$train$paths$n_paths())
  pes1 <- event_risk(p_mfit, "Death")

  # Test that solving time evolution with single draw works
  tp2 <- solve_trans_prob_fit(mfit)
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
  expect_true(inherits(fit, "MultistateModelFit"))

  # Plot baseline hazards
  p0 <- fit$plot_h0()
  expect_s3_class(p0, "ggplot")

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
  p_oos <- generate_paths(fit, data = jd$test)
  P <- solve_trans_prob_fit(fit)
  P_oos <- solve_trans_prob_fit(fit, data = jd$test)
  expect_equal(nrow(P), 75)
  expect_equal(nrow(P_oos), 25)
})
