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
  expect_s3_class(fit$plot_basisfun(), "ggplot")
  expect_s3_class(fit$plot_h0(), "ggplot")

  # Pathgen
  p <- generate_paths(fit, n_rep = 3)
  expect_true(inherits(p, "PathData"))

  # P(event)
  pe <- p_event(p)
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
  expect_s3_class(fit$plot_h0(), "ggplot")

  # Pathgen
  NR <- 4
  p <- generate_paths(fit, n_rep = NR)
  expect_true(inherits(p, "PathData"))
  expect_equal(p$n_paths(), options$iter_sampling * options$N_subject * NR)

  # P(event)
  pe <- p_event(p)
  expect_equal(nrow(pe), 2)

  # Solve
  P <- solve_trans_prob_matrix(mod$system, c(-7, -7, -7))
  checkmate::assert_matrix(P, ncols = 3, nrows = 3)
  r <- solve_trans_prob_fit(fit)
  checkmate::assert_data_frame(r, ncols = 4)
})
