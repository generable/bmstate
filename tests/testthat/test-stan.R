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
  jd <- mod$simulate_data(options$N_subject, w0 = exp(-4))

  # Fit
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains
  )

  expect_true(inherits(fit, "MultistateModelStanFit"))

  # Plots
  expect_s3_class(fit$plot_basisfun(), "ggplot")
  expect_s3_class(fit$plot_h0(), "ggplot")

  # Pathgen
  p <- generate_paths(fit)
  expect_true(inherits(p, "PathData"))
})

test_that("fitting with Stan works (multi-transition)", {
  # Setup
  tm <- transmat_illnessdeath()
  mod <- create_msm(tm)
  jd <- mod$simulate_data(options$N_subject, w0 = exp(-3))

  # Fit
  fit <- fit_stan(mod, jd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains
  )

  expect_true(inherits(fit, "MultistateModelStanFit"))
  expect_s3_class(fit$plot_h0(), "ggplot")

  # Pathgen
  p <- generate_paths(fit)
  expect_true(inherits(p, "PathData"))
})
