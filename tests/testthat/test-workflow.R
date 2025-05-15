test_that("entire workflow works", {
  # Options
  options <- list(
    N_subject = 200,
    iter_warmup = 60,
    iter_sampling = 30,
    chains = 1,
    NK = 3
  )
  NP <- options$iter_sampling
  NR <- 2

  # Setup
  h0_base <- 1e-3
  setup <- example_sim_setup_illnessdeath()

  # Simulate data
  pd <- setup$mod$simulate_data(options$N_subject, beta_haz = setup$beta_haz)

  split <- do_split(pd)
  expect_true(inherits(pd, "PathData"))
  expect_equal(pd$n_paths(), options$N_subject)
  expect_equal(pd$longest_path()$n_paths(), 1)

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
  sub_ids_char <- pd$subject_df |> dplyr::filter(subject_id %in% split$train_sub)
  sub_ids_char_test <- pd$subject_df |> dplyr::filter(subject_id %in% split$test_sub)

  # Paths from 0 to 3 years starting from Randomization
  paths_3yr <- generate_paths_many_subjects(gq, sd, possible_covs,
    TFI, init_state, t_start, t_max_gen,
    state_names = pd$state_names,
    trans_names = dt$legend$trans_char,
    df_subjects = sub_first_rows,
    terminal_states = pd$terminal_states,
    num_paths = NP, n_repeats = NR,
  )
  expect_equal(paths_3yr$n_paths(), 9000)

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
  expect_equal(paths_3yr_oos$n_paths(), 3000)

  h0_true_all <- d$h0
  m_sub <- d$m_sub
  df_beta_true <- d$df_beta_true
  N_subject <- options$N_subject
  res <- dplyr::lst(
    fit, stan_dat, pd, paths_3yr, paths_3yr_oos,
    sub_first_rows, sub_first_rows_test, N_subject, df_beta_true,
    h0_true_all, m_sub, split
  )

  # Covariate effect and baseline hazard plots
  sd <- res$stan_dat$stan_data
  all_states <- res$pd$state_names
  dt <- res$pd$as_transitions()
  plt_oth <- plot_other_beta(res$fit, res$stan_dat, res$pd, res$df_beta_true)
  expect_equal(length(plt_oth), 3)

  df_h0 <- tibble::as_tibble(data.frame(
    transition = seq_len(length(res$h0_true_all)),
    log_h0_true = log(res$h0_true_all)
  ))
  plt_h0 <- plot_h0(res$fit, sd, res$pd, df_h0)
  expect_s3_class(plt_h0, "ggplot")

  # Event probabilities
  names <- c("1. Obs", "2. IS pred", "3. OOS pred")
  plt_a <- plot_p_event(res$pd, res$paths_3yr, res$paths_3yr_oos, names)

  # no points selected for one or more curves, consider using the extend argument
  plt_b <- plot_p_event_by_dose(res$pd, res$paths_3yr, res$paths_3yr_oos, names)
  expect_s3_class(plt_a, "ggplot")
  expect_s3_class(plt_b, "ggplot")

  # Event summary
  # There was 1 warning in `dplyr::mutate()`.
  # i In argument: `event_surv = furrr::future_map(...)`.
  # Caused by warning in `serializedSize()`:
  #  ! 'package:stats' may not be available when loading
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
  expect_true(inherits(ppsurv_subj, "tbl_df"))
  expect_true(inherits(ppsurv_subj_oos, "tbl_df"))

  # Brier score
  train_sub <- res$sub_first_rows$subject_id
  test_sub <- res$sub_first_rows_test$subject_id
  plt_bs_is <- create_brier_score_plot(
    ppsurv_subj, res$pd, train_sub
  )
  plt_bs_oos <- create_brier_score_plot(
    ppsurv_subj_oos, res$pd, test_sub
  )
  expect_s3_class(plt_bs_is, "ggplot")
  expect_s3_class(plt_bs_oos, "ggplot")

  # Concordance index
  ci_eval <- create_cindex_plot(
    res$pd, res$paths_3yr, res$paths_3yr_oos,
    train_sub, test_sub,
    ppsurv_subj, ppsurv_subj_oos
  )
  expect_equal(length(ci_eval), 4)
  expect_s3_class(ci_eval$plt_ci_is, "ggplot")
  expect_s3_class(ci_eval$plt_ci_oos, "ggplot")
})
