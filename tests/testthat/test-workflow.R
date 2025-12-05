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

  # Test
  a <- as_any_event(jd$train$paths, "Healthy")
  checkmate::assert_true(all(a$state_names() == c("Healthy", "Any event")))

  # Print
  expect_output(print(fit), "A MultistateModelFit with 30 draws")

  # Debug functions
  a <- plot_stan_data_integral(mod, fit$get_data(), 1)
  expect_true(is_ggplot(a))
  a <- plot_stan_data_matrix(mod, fit$get_data(), "transition", 1)
  expect_true(is_ggplot(a))
  expect_error(plot_stan_data_matrix(mod, fit$get_data(), "asfd", 4))

  # Covariate effects
  beta <- fit$covariate_effects()
  expect_equal(nrow(beta), 6)

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
    refresh = 20,
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
  # Create model
  sn <- c("Alive", "Death")
  tm <- transmat_survival(state_names = sn)
  t3yr <- 3 * 365.25
  covs <- c("age", "sex")
  pk_covs <- list(
    CL = "CrCL",
    V2 = c("weight", "sex")
  )
  mod <- create_msm(tm, covs, pk_covs,
    num_knots = 5, t_max = t3yr,
    categ_covs = "sex"
  )
  bh_true <- matrix(0, 1, 3)
  bh_true[1, 1] <- 0.5 # age effect on death
  bh_true[1, 2] <- -0.5 # xpsr effect on death
  sn <- "Death"
  rownames(bh_true) <- paste0("Effect on ", sn)
  colnames(bh_true) <- mod$covs()

  h0_true <- 1e-3
  beta_pk <- list(CL = -0.3, V2 = c(0.2, 0.3))

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
    refresh = 20,
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
  pk_covs <- list(
    V2 = "weight"
  )
  tm <- transmat_survival()
  tmax <- 3 * 365.25
  mod <- create_msm(tm, hazard_covs = NULL, pk_covs, t_max = tmax, pk_only = TRUE)
  jd <- mod$simulate_data(options$N_subject)

  # Fit the model
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains,
    refresh = 20,
    adapt_delta = 0.95,
    init = 0.1
  )

  # Plot pk fit
  sid <- unique(jd$paths$subject_df$subject_id)[1:6]
  plt <- fit$plot_pk(subject_ids = sid)
  ttt <- seq(100, 300, by = 1)
  tl <- list()
  for (j in 1:options$N_subject) {
    tl[[j]] <- seq(180, 220, by = 1)
  }
  dd <- fit$data$dosing
  TH_true <- as.matrix(fit$data$paths$subject_df[, c("ka", "CL", "V2")])
  TH_fit <- as.matrix(msmfit_pk_params(fit$mean_fit())[[1]])
  pks1 <- dd$simulate_pk(tl, TH_true, 10^5)
  pks1 <- pks1 |> dplyr::filter(.data$subject_id %in% sid)
  pks2 <- dd$simulate_pk(tl, TH_fit, 10^5)
  pks2 <- pks2 |> dplyr::filter(.data$subject_id %in% sid)
  plt2 <- plt + geom_line(data = pks1, aes(x = time, y = val, group = subject_id)) +
    geom_line(data = pks2, aes(x = time, y = val, group = subject_id), color = "red")
  expect_true(is_ggplot(plt2))

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
