test_that("cav data analysis works", {
  library(ggplot2)
  library(dplyr)

  # Modify original data like in
  # https://rviews.rstudio.com/2023/04/19/multistate-models-for-medical-applications/
  df_full <- as_tibble(msm::cav) |>
    group_by(PTNUM) |>
    mutate(
      min_age = min(age),
      state = cummax(state)
    )
  df_full$PTNUM <- as.character(df_full$PTNUM)

  # Simplify the snapshots data more
  df_snapshots <- df_full |> select(PTNUM, state, years, min_age, dage, sex)
  df_snapshots <- df_snapshots |> rename(
    subject_id = PTNUM,
    age = min_age,
    time = years
  )

  # Collapse
  df_state_changes <- df_snapshots |>
    arrange(subject_id, time) |>
    group_by(subject_id) |>
    mutate(
      .first = row_number() == 1,
      .last  = row_number() == n(),
      .chg   = state != lag(state) # TRUE at the first row of each new state run
    ) |>
    filter(.first | .last | .chg | is.na(lag(state))) |> # keep first, last, and changes
    select(-.first, -.last, -.chg) |>
    ungroup()

  # Every row is a transition but the first
  df_state_changes <- df_state_changes |>
    group_by(subject_id) |>
    mutate(is_transition = tidyr::replace_na(state != lag(state), FALSE)) |>
    ungroup()

  # To PathData
  tm <- transmat_progression()
  covs <- c("age")
  pd <- df_to_pathdata(df_state_changes, tm, covs)
  pd0 <- as_single_event(pd, "Dead")
  tm0 <- pd0$transmat
  t_ev <- pd0$transition_times()

  # Create model and fit
  mod0 <- create_msm(tm0,
    hazard_covs = covs,
    n_grid = 10000
  )
  mod0$set_knots(t_max = 20, t_event = t_ev, num_knots = 3)
  fit0 <- fit_stan(mod0,
    data = pd0,
    chains = 1,
    iter_sampling = 30,
    iter_warmup = 30
  )

  # Plot baseline hazard
  plt0 <- fit0$plot_h0()
  expect_true(is_ggplot(plt0))

  # State occupancy
  f0 <- fit0$mean_fit()
  a <- p_state_occupancy(f0, data = f0$data, oos = TRUE)
  plt <- plot_state_occupancy(a)
  expect_true(is_ggplot(plt))
})
