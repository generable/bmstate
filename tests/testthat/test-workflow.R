test_that("entire workflow works", {
  # Options
  options <- list(
    N_subject = 100,
    covariates = c("sex", "first_dose_amount", "age"),
    iter_warmup = 60,
    iter_sampling = 30,
    chains = 1,
    NK = 3
  )
  NP <- options$iter_sampling
  NR <- 2

  # Setup
  effect_sizes <- c(0.5, 0.5, 0)
  h0_base <- 1e-3

  # Create model
  mod <- create_stan_model()
  expect_true(inherits(mod, "CmdStanModel"))

  # Simulate data
  d <- simulate_example_data(options$N_subject,
    sys_idx = 1, h0_base = h0_base,
    effect_sizes = effect_sizes
  )
  possible_covs <- d$covs
  pd <- d$pd
  split <- do_split(pd)
  expect_true(inherits(pd, "PathData"))
  expect_equal(pd$n_paths(), options$N_subject)

  # Stan data
  stan_dat <- create_stan_data(
    pd, NULL, options$covariates, "age", "age", "age",
    split$train_sub, split$test_sub,
    P = 30, NK = options$NK,
    do_pk = FALSE,
    prior_only = FALSE
  )
  sd <- stan_dat$stan_data

  # Fit the model
  fit <- mod$sample(
    data = sd,
    iter_warmup = options$iter_warmup,
    iter_sampling = options$iter_sampling,
    chains = options$chains,
    refresh = 500,
    adapt_delta = 0.95,
    init = 0.1
  )

  # Path prediction
  gq <- mod$generate_quantities(fit, data = sd)
  dt <- pd$as_transitions()
  TFI <- legend_to_TFI_matrix(dt$legend)
  t_max_gen <- 3 * 365.25 # 3 yr
  t_start <- 0
  init_state <- 1

  # Finding subjects for whom to predict
  sub_first_rows <- subject_df_with_idx(pd, split$train_sub, stan_dat$id_map_train)
  sub_first_rows_test <- subject_df_with_idx(pd, split$test_sub, stan_dat$id_map_test)
  sub_ids_char <- pd$subject_df |> filter(subject_id %in% split$train_sub)
  sub_ids_char_test <- pd$subject_df |> filter(subject_id %in% split$test_sub)

  # Paths from 0 to 3 years starting from Randomization
  paths_3yr <- generate_paths_many_subjects(gq, sd, possible_covs,
    TFI, init_state, t_start, t_max_gen,
    state_names = pd$state_names,
    trans_names = dt$legend$trans_char,
    df_subjects = sub_first_rows,
    terminal_states = pd$terminal_states,
    num_paths = NP, n_repeats = NR,
  )
  expect_equal(paths_3yr$n_paths(), 900)

  # Paths from 0 to 3 years starting from Randomization
  paths_3yr_oos <- generate_paths_many_subjects(gq, sd, possible_covs,
    TFI, init_state, t_start, t_max_gen,
    state_names = pd$state_names,
    trans_names = dt$legend$trans_char,
    df_subjects = sub_first_rows_test,
    terminal_states = pd$terminal_states,
    num_paths = NP, n_repeats = NR,
    oos = TRUE
  )
  expect_equal(paths_3yr_oos$n_paths(), 300)

  h0_true_all <- d$h0
  m_sub <- d$m_sub
  df_beta_true <- d$df_beta_true
  res <- dplyr::lst(
    fit, stan_dat, pd, paths_3yr, paths_3yr_oos,
    sub_first_rows, sub_first_rows_test, N_subject, df_beta_true,
    h0_true_all, m_sub, split
  )

  # Study the covariate effects and baseline hazard
  sd <- res$stan_dat$stan_data
  all_states <- res$pd$state_names
  dt <- res$pd$as_transitions()
  plt_oth <- NULL
  if (sd$N_oth > 0) {
    plt_oth <- plot_other_beta(res$fit, res$stan_dat, res$pd, res$df_beta_true)
    ff <- function(x) {
      x + theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
    }
    plt_oth <- lapply(plt_oth, ff)
  }
  df_h0 <- tibble::as_tibble(data.frame(
    transition = seq_len(length(res$h0_true_all)),
    log_h0_true = log(res$h0_true_all)
  ))
  ttl <- paste0("time = ", round(res$fit$time()$total, 3))
  plt_h0 <- plot_h0(res$fit, sd, res$pd, df_h0) +
    xlab("Time (days)") + theme(legend.position = "none") +
    ggtitle(ttl)
  list_beta <- c(plt_oth)
  plt_beta <- ggarrange(plotlist = list_beta)

  # Event probabilities
  names <- c("1. Obs", "2. IS pred", "3. OOS pred")
  plt_a <- plot_p_event(res$pd, res$paths_3yr, res$paths_3yr_oos, names)
  plt_b <- plot_p_event_by_dose(res$pd, res$paths_3yr, res$paths_3yr_oos, names) +
    ggtitle("By dose")
  plt_freq <- ggarrange(plt_a, plt_b, nrow = 1, ncol = 2)

  # Brier score
  target_times <- seq(0, 3 * 365.25, by = 365.25 / 4)
  ppsurv_subj <- summarize_ppsurv(
    res$paths_3yr,
    target_times = target_times,
    by = c("subject_id"),
    truncate_at_terminal_events = FALSE
  )
  ppsurv_subj_oos <- summarize_ppsurv(
    res$paths_3yr_oos,
    target_times = target_times,
    by = c("subject_id"),
    truncate_at_terminal_events = FALSE
  )
  train_sub <- res$sub_first_rows$subject_id
  test_sub <- res$sub_first_rows_test$subject_id
  plt_bs_is <- create_brier_score_plot(
    ppsurv_subj, res$pd, train_sub
  ) + ggtitle("Brier scores (in-sample)")
  plt_bs_oos <- create_brier_score_plot(
    ppsurv_subj_oos, res$pd, test_sub
  ) + ggtitle("Brier scores (out-of-sample)")
  plt_brier_all <- ggarrange(plt_bs_is, plt_bs_oos, nrow = 2, ncol = 1)

  # Concordance index
  ci_eval <- create_cindex_plot(
    res$pd, res$paths_3yr, res$paths_3yr_oos,
    train_sub, test_sub,
    ppsurv_subj, ppsurv_subj_oos
  )
  plt_ci_all <- ggarrange(ci_eval$plt_ci_is, ci_eval$plt_ci_oos, nrow = 2, ncol = 1)
})
