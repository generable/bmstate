# Options
options <- list(
  N_subject = 100,
  iter_warmup = 100,
  iter_sampling = 100,
  chains = 1
)

test_that("fitting with Stan works (single transition)", {
  # Setup
  h0_base <- 1e-3
  tm <- transmat_survival()
  mod <- create_msm(tm)

  # Simulate data
  jd <- mod$simulate_data(options$N_subject)

  # Fit
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains
  )

  expect_true(inherits(fit, "MultistateModelStanFit"))
})

test_that("fitting with Stan works (multi-transition)", {
  # Setup
  tm <- transmat_illnessdeath()
  mod <- create_msm(tm)
  jd <- mod$simulate_data(options$N_subject, w0 = 1e-2)

  # Fit
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains
  )

  expect_true(inherits(fit, "MultistateModelStanFit"))
})
