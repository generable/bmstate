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

  # Test covariates
  a <- potential_covariates(jd$train$paths)
  expect_equal(colnames(a), c("pval", "target_state", "covariate"))

  # CoxPH fit
  cph <- jd$train$paths$fit_coxph(covariates = mod$covs())
  msf <- jd$train$paths$fit_mstate()
  expect_true(inherits(msf, "msfit"))

  # Fit the model
  a <- fit_stan(mod, jd$train,
    refresh = 5,
    init = 0.1,
    return_stanfit = TRUE,
    method = "pathfinder",
    draws = options$iter_sampling
  )
  fit <- a$fit
  expect_true(inherits(fit, "MultistateModelFit"))

  # Computing hazard multipliers
  log_m <- msmfit_log_hazard_multipliers(fit)
  log_m_test <- msmfit_log_hazard_multipliers(fit, oos = TRUE, data = jd$test)

  # Path generation
  p <- generate_paths(fit, n_rep = 3)
  p_oos <- generate_paths(fit, n_rep = 3, oos = TRUE, data = jd$test)

  # Path generation starting from later time
  p1 <- generate_paths(fit, n_rep = 3, t_start = 600)

  pe <- p_state_visit(p)
  pe1 <- p_state_visit(p1)
  pe_oos <- p_state_visit(p_oos)
  pe_bysex <- p_state_visit(p, by = "sex")
  pe_bydose <- p_state_visit(p, by = "dose_amt")
  expect_equal(nrow(pe), 2)
  expect_equal(nrow(pe1), 2)
  expect_equal(nrow(pe_bysex), 4)
  expect_equal(nrow(pe_bydose), 6)

  # Baseline hazard viz
  plot_h0 <- fit$plot_h0() # should be at h0_base level
  expect_true(is_ggplot(plot_h0))

  # Scoring
  ev <- "Dead"
  a <- create_scoring_df(jd$train$paths, p, ev)

  ci <- c_index(a)
  expect_gt(ci$concordance, 0)

  # Test that reducing to mean works
  mfit <- fit$mean_fit()
  expect_equal(mfit$num_draws(), 1)
  p_mfit <- generate_paths(mfit, n_rep = 100)
  expect_equal(p_mfit$n_paths(), 100 * jd$train$paths$n_paths())
  pes1 <- p_state_visit_per_subject(p_mfit, ev)

  # Test that solving time evolution with single draw works
  tp2 <- p_state_occupancy(mfit)
  pp <- plot_state_occupancy(tp2)
  expect_true(is_ggplot(pp))

  # Single-transition model
  tds <- as_single_event(jd$train$paths, event = ev)
  tds <- JointData$new(tds, jd$train$dosing)
  mod_tte <- create_msm(tds$paths$transmat)
  a <- fit_stan(mod_tte, tds,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains,
    refresh = 5,
    adapt_delta = 0.95,
    init = 0.1,
    return_stanfit = TRUE
  )
  fit_tte <- a$fit$mean_fit()
  r <- p_state_occupancy(fit_tte)
  pp <- plot_state_occupancy(r)
  expect_true(is_ggplot(pp))
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
  fit <- fit$mean_fit()
  expect_true(inherits(fit, "MultistateModelFit"))

  # Plot baseline hazards
  p0 <- fit$plot_h0()
  expect_true(is_ggplot(p0))

  # PK params
  pkpar <- msmfit_pk_params(fit)
  pkpar_oos <- msmfit_pk_params(fit, oos = TRUE, data = jd$test)
  expect_true(length(pkpar) == 1)
  expect_true(length(pkpar_oos) == 1)

  # Hazard multipliers
  log_m <- msmfit_log_hazard_multipliers(fit)
  log_m_test <- msmfit_log_hazard_multipliers(fit, oos = TRUE, jd$test)
  N_sub_test <- length(jd$test$paths$unique_subjects())
  H <- mod$system$num_trans()
  expect_equal(dim(log_m_test[[1]]), c(N_sub_test, H))

  # Path simulation
  p <- generate_paths(fit)
  p_oos <- generate_paths(fit, oos = TRUE, data = jd$test)
  P <- p_state_occupancy(fit)
  P_oos <- p_state_occupancy(fit, oos = TRUE, data = jd$test)
  expect_equal(nrow(P), 75 * 2 * 30)
  expect_equal(nrow(P_oos), 25 * 2 * 30)

  # PK fit plot
  pf1 <- fit$plot_pk()
  pf2 <- fit$plot_pk(oos = TRUE, data = jd$test)
  expect_true(is_ggplot(pf1))
  expect_true(is_ggplot(pf2))
})

test_that("entire workflow works (with single transition PK)", {
  # Setup
  h0_base <- 1e-3

  # Create model
  sn <- c("Alive", "Death")
  tm <- transmat_survival(state_names = sn)
  t3yr <- 3 * 365.25
  covs <- c("age")
  pk_covs <- list(
    CL = "CrCL",
    V2 = c("weight", "sex")
  )
  mod <- create_msm(tm, covs, pk_covs, num_knots = 5, t_max = t3yr)
  bh_true <- matrix(0, 1, 2)
  bh_true[1, 1] <- 0.5 # age effect on death
  bh_true[1, 2] <- -0.5 # dose effect on death
  sn <- "Death"
  rownames(bh_true) <- paste0("Effect on ", sn)
  colnames(bh_true) <- mod$covs()

  h0_true <- 1e-3
  beta_pk <- list(CL = -0.1, V2 = c(0, 0.1))

  # Simulate data
  simdat <- mod$simulate_data(
    options$N_subject,
    beta_haz = bh_true,
    beta_pk = beta_pk,
    w0 = h0_true
  )

  # Split
  jd <- simdat
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
  fit <- fit$mean_fit()
  expect_true(inherits(fit, "MultistateModelFit"))

  # Plot baseline hazards
  p0 <- fit$plot_h0()
  expect_true(is_ggplot(p0))

  # PK params
  pkpar <- msmfit_pk_params(fit)
  pkpar_oos <- msmfit_pk_params(fit, oos = TRUE, jd$test)
  expect_true(length(pkpar) == 1)
  expect_true(length(pkpar_oos) == 1)

  # Hazard multipliers
  log_m <- msmfit_log_hazard_multipliers(fit)
  log_m_test <- msmfit_log_hazard_multipliers(fit, oos = TRUE, data = jd$test)
  N_sub_test <- length(jd$test$paths$unique_subjects())
  H <- mod$system$num_trans()
  expect_equal(dim(log_m_test[[1]]), c(N_sub_test, H))
})


test_that("PK-only works", {
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
    init = 0.1,
    pk_only = TRUE
  )
  fit <- fit$mean_fit()
  expect_true(inherits(fit, "MultistateModelFit"))

  # Plot baseline hazards
  expect_error(fit$plot_h0(), "This is a PK-only fit")

  # PK params
  pkpar <- msmfit_pk_params(fit)
  pkpar_oos <- msmfit_pk_params(fit, oos = TRUE, data = jd$test)
  expect_true(length(pkpar) == 1)
  expect_true(length(pkpar_oos) == 1)

  # PK fit plot
  pf1 <- fit$plot_pk()
  pf2 <- fit$plot_pk(data = jd$test, oos = TRUE)
  expect_true(is_ggplot(pf1))
  expect_true(is_ggplot(pf2))
})
