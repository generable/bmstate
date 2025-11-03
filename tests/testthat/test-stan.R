# Options
options <- list(
  N_subject = 30,
  iter_warmup = 100,
  iter_sampling = 20,
  chains = 1
)

test_that("fitting with Stan works (single transition)", {
  # Setup
  tm <- transmat_survival()
  mod <- create_msm(tm)

  # Simulate data
  jd <- mod$simulate_data(options$N_subject, w0 = exp(-5))

  # Fit
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains
  )

  expect_true(inherits(fit, "MultistateModelFit"))

  # Plots
  expect_true(is_ggplot(fit$plot_basisfun()))
  expect_true(is_ggplot(fit$plot_h0()))

  # Pathgen
  p <- generate_paths(fit, n_rep = 3)
  expect_true(inherits(p, "PathData"))

  # P(event)
  pe <- p_state_visit(p)
  expect_equal(nrow(pe), 1)
})

test_that("fitting with Stan works (multi-transition)", {
  # Setup
  tm <- transmat_illnessdeath()
  mod <- create_msm(tm)
  jd <- mod$simulate_data(options$N_subject, w0 = exp(-6))

  # Fit
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains
  )

  expect_true(inherits(fit, "MultistateModelFit"))
  expect_true(is_ggplot(fit$plot_h0()))

  # Pathgen
  NR <- 4
  p <- generate_paths(fit, n_rep = NR)
  expect_true(inherits(p, "PathData"))
  expect_equal(p$n_paths(), options$iter_sampling * options$N_subject * NR)

  # P(event)
  pe <- p_state_visit(p)
  expect_equal(nrow(pe), 2)

  # Solve
  P <- solve_trans_prob_matrix(mod$system, c(0, 1), c(-7, -7, -7))
  expect_equal(dim(P), c(2, 3, 3))
  r <- p_state_occupancy(fit)
  pp <- plot_state_occupancy(r)
  expect_true(is_ggplot(pp))

  checkmate::assert_data_frame(r, ncols = 4)
})
