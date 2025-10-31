library(msm)
library(dplyr)

test_that("cav data analysis works", {
  # Modify original data like in
  # https://rviews.rstudio.com/2023/04/19/multistate-models-for-medical-applications/
  df_full <- as_tibble(cav) |>
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

  # To PathData
  tm <- transmat_progression()
  covs <- c("age", "dage", "sex")
  pd <- df_to_pathdata(df_state_changes, tm, covs)
  pd0 <- as_single_event(pd, "Dead")
  tm0 <- pd0$transmat

  mod0 <- create_msm(tm0,
    hazard_covs = covs,
    categ_covs = "sex",
    num_knots = 5,
    tmax = 20,
    n_grid = 1000
  )
  dat0 <- JointData$new(pd0, NULL)
  fit0 <- fit_stan(mod0,
    data = dat0, chains = 1, iter_sampling = 100,
    iter_warmup = 100
  )
  fit0_pm <- fit0$mean_fit()

  # Plot baseline hazard
  plt0 <- fit0$plot_h0()
})
