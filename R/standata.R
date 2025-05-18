#' Creating Stan data list
#'
#' @export
#' @inheritParams fit_model
#' @return A list of data for Stan.
create_stan_data <- function(model, pd, dosing = NULL, prior_only = FALSE,
                             delta_grid = 1) {
  checkmate::assert_class(model, "MultistateModel")
  checkmate::assert_class(pd, "PathData")
  checkmate::assert_number(delta_grid, lower = 0)
  tm <- pd$transmat
  check_equal_transmats(tm, model$system$tm())
  checkmate::assert_logical(prior_only, len = 1)
  if (!is.null(dosing)) {
    checkmate::assert_class(dosing, "data.frame")
    stopifnot(!all(duplicated(dosing$subject_id)))
  }

  # Initial Stan data
  stan_dat <- c(
    create_stan_data_idx_sub(pd),
    create_stan_data_indicators(pd),
    create_stan_data_spline(pd, model, delta_grid),
    create_stan_data_covariates(pd, model),
    create_stan_data_trans_types(pd)
  )

  # PK data
  out <- c(out, create_stan_data_pk(pd, model))
  out <- c(out, create_stan_data_time_since_last_pk(out))

  # Likelihood flags
  out$omit_lik_hazard <- as.integer(prior_only)
  out$omit_lik_pk <- as.integer(prior_only)
  out$do_pk <- as.integer(do_pk)

  # Add t_int if needed by some model versions
  out$t_int <- out$t_end - out$t_start

  # Return
  list(
    stan_data = out,
    id_map_train = id_map_train,
    id_map_test = id_map_test,
    x_oth_names = oth$x_oth_names
  )
}

# Transition types
create_stan_data_trans_types <- function(pd) {
  tm <- pd$transmat
  # Average event rates
  df_ttype <- average_haz_per_ttype(pd) |> dplyr::arrange(.data$trans_idx)
  ttype <- tm$trans_df()$trans_type
  list(
    N_trans_types = max(ttype),
    ttype = ttype,
    mu_w0 = df_ttype$log_h0_avg
  )
}

# Create binary indicator matrices
create_stan_data_indicators <- function(pd) {
  tm <- pd$transmat
  dat <- pd$as_transitions()
  N_int <- nrow(dat)
  N_trans <- tm$num_trans()
  transition <- matrix(0, N_trans, N_int)
  at_risk <- matrix(0, N_trans, N_int)
  for (n in seq_len(N_int)) {
    tval <- dat$trans_idx[n]
    if (tval > 0) {
      transition[tval, n] <- 1
    }
    at_risk[tm$at_risk(dat$from[n]), n] <- 1
  }
  out <- list(
    at_risk = at_risk,
    transition = transition,
    N_int = N_int,
    N_trans = N_trans
  )

  # Another way to format the transitions and at risk matrices
  # needed for vectorizing likelihood
  out <- c(out, which_format_for_stan(transition, "trans"))
  out <- c(out, which_format_for_stan(at_risk, "risk"))

  # Return
  out
}

# Which subject does each interval correspond to
create_stan_data_idx_sub <- function(pd) {
  idx_sub <- as.integer(pd$as_transitions()$subject_index)
  N_sub <- max(idx_sub)
  stopifnot(all(seq_len(N_sub) %in% idx_sub))
  list(
    N_sub = N_sub,
    idx_sub = idx_sub
  )
}

# Evaluate spline basis functions for each interval
create_stan_data_spline <- function(pd, model, delta_grid) {
  # Interval end time spline evaluations
  dat <- pd$as_transitions()
  t <- dat$time
  SBF <- model$system$basisfun_matrix(t)
  t_max <- model$system$get_tmax()
  if (delta_grid > 10 * t_max) {
    stop("delta_grid is very large compared to t_max")
  }

  # Grid spline evaluations
  t_grid <- seq(delta_grid / 2, ceiling(t_max), by = delta_grid) # midpoints
  SBF_grid <- model$system$basisfun_matrix(t_grid)

  # Grid interval indices
  sd_int <- create_stan_data_intervalidx(
    dat$time_prev, dat$time, t_grid, delta_grid
  )

  # Return
  sd <- list(
    N_sbf = ncol(SBF),
    SBF = SBF,
    SBF_grid = SBF_grid,
    N_grid = length(t_grid),
    t_grid = t_grid,
    delta_grid = delta_grid
  )
  c(sd, sd_int)
}

# Get original subject id based on numeric id that is in Stan data
subject_idx_to_id <- function(id_map, idx) {
  id_map$subject_id[which(id_map$x_sub == idx)]
}

# PK data
create_stan_data_pk <- function(pk, sd, id_map_train, id_map_test) {
  N_sub <- sd$N_sub
  N_sub_oos <- sd$N_sub_oos
  t_obs_pk <- matrix(0, N_sub, 2)
  t_obs_pk_oos <- matrix(0, N_sub_oos, 2)
  conc_pk <- matrix(0, N_sub, 2)
  conc_pk_oos <- matrix(0, N_sub_oos, 2)
  if (is.null(pk)) {
    last_n <- 2
  } else {
    last_n <- unique(pk$pk_last_n)
  }
  checkmate::assert_integerish(last_n, lower = 2, len = 1)
  last_times <- matrix(0, N_sub, last_n)
  last_doses <- matrix(0, N_sub, last_n)
  last_two_times <- matrix(0, N_sub, 2)
  last_two_doses <- matrix(0, N_sub, 2)
  last_times_oos <- matrix(0, N_sub_oos, last_n)
  last_doses_oos <- matrix(0, N_sub_oos, last_n)
  pk_lloq <- rep(0, N_sub)
  for (n in 1:N_sub) {
    sub_id <- subject_idx_to_id(id_map_train, n)
    if (is.null(pk)) {
      conc_pk[n, ] <- 1
      t_obs_pk[n, ] <- 1
      last_times[n, ] <- 1
      last_doses[n, ] <- 1
      last_two_times[n, ] <- 1
      last_two_doses[n, ] <- 1
      pk_lloq[n] <- 1
    } else {
      pk_sub <- pk |> dplyr::filter(subject_id == sub_id)
      if (nrow(pk_sub) != 1) {
        stop("found multiple rows for subject ", n)
      }
      conc_pk[n, ] <- unlist(pk_sub$list_pk_conc)
      t_obs_pk[n, ] <- unlist(pk_sub$list_pk_hours)
      last_times[n, ] <- unlist(pk_sub$list_dose_hours)
      last_doses[n, ] <- unlist(pk_sub$list_dose_amts)
      last_two_times[n, ] <- last_times[n, (last_n - 1):last_n]
      last_two_doses[n, ] <- last_doses[n, (last_n - 1):last_n]
      pk_lloq[n] <- pk_sub$pk_lloq
    }
  }
  for (n in 1:N_sub_oos) {
    sub_id <- subject_idx_to_id(id_map_test, n)
    if (is.null(pk)) {
      conc_pk_oos[n, ] <- 1
      t_obs_pk_oos[n, ] <- 1
      last_times_oos[n, ] <- 1
      last_doses_oos[n, ] <- 1
    } else {
      pk_sub <- pk |> dplyr::filter(subject_id == sub_id)
      if (nrow(pk_sub) != 1) {
        stop("found multiple rows for subject ", n)
      }
      conc_pk_oos[n, ] <- unlist(pk_sub$list_pk_conc)
      t_obs_pk_oos[n, ] <- unlist(pk_sub$list_pk_hours)
      last_times_oos[n, ] <- unlist(pk_sub$list_dose_hours)
      last_doses_oos[n, ] <- unlist(pk_sub$list_dose_amts)
    }
  }

  # Return
  list(
    t_obs_pk = t_obs_pk,
    conc_pk = conc_pk,
    t_obs_pk_oos = t_obs_pk_oos,
    conc_pk_oos = conc_pk_oos,
    dose_ss = sd$x_fda[sd$sub_start_idx],
    dose_ss_oos = sd$x_fda_oos[sd$sub_start_idx_oos],
    t_pred_pk = t_obs_pk,
    t_pred_pk_oos = t_obs_pk_oos,
    pred_pk = 0, # will be set to 1 when using standalone GQ
    N_pred_pk = 2,
    N_last = last_n,
    last_times = last_times,
    last_doses = last_doses,
    last_two_times = last_two_times,
    last_two_doses = last_two_doses,
    last_times_oos = last_times_oos,
    last_doses_oos = last_doses_oos,
    pk_lloq = pk_lloq,
    I_auc = as.numeric(model$has_pk())
  )
}

# PK data
update_stan_data_pk_pred <- function(sd, P = 200) {
  N_sub <- sd$N_sub
  N_sub_oos <- sd$N_sub_oos
  t_pred_pk <- matrix(0, N_sub, P)
  t_pred_pk_oos <- matrix(0, N_sub_oos, P)
  for (n in 1:N_sub) {
    t_max <- 1.002 * (max(sd$last_times[n, ]) + 1 * 24)
    t_min <- 0.998 * (min(sd$last_times[n, ]) - 2 * 24)
    t_pred_pk[n, ] <- seq(t_min, t_max, length.out = P)
  }
  for (n in 1:N_sub_oos) {
    t_max <- 1.002 * (max(sd$last_times_oos[n, ]) + 1 * 24)
    t_min <- 0.998 * (min(sd$last_times_oos[n, ]) - 2 * 24)
    t_pred_pk_oos[n, ] <- seq(t_min, t_max, length.out = P)
  }
  sd$t_pred_pk <- t_pred_pk
  sd$t_pred_pk_oos <- t_pred_pk_oos
  sd$pred_pk <- 1
  sd$N_pred_pk <- P
  sd
}

# Time since last dose
create_stan_data_time_since_last_pk <- function(sd) {
  N_sub <- sd$N_sub
  N_sub_oos <- sd$N_sub_oos
  t_since_last_pk <- matrix(0, N_sub, 2)
  t_since_last_pk_oos <- matrix(0, N_sub_oos, 2)
  for (n in 1:N_sub) {
    t_since_last_pk[n, ] <- time_since_last_dose(
      sd$t_obs_pk[n, ], sd$last_times[n, ]
    )
  }
  for (n in 1:N_sub_oos) {
    t_since_last_pk_oos[n, ] <- time_since_last_dose(
      sd$t_obs_pk_oos[n, ], sd$last_times_oos[n, ]
    )
  }
  list(
    t_since_last_pk = t_since_last_pk,
    t_since_last_pk_oos = t_since_last_pk_oos
  )
}

# Normalized covariates to stan data
standata_scaled_covariates <- function(pd, model, name) {
  covs <- model$data_covs()
  sub_df <- pd$subject_df[, c("subject_index", covs)]
  x <- list()
  nc <- length(covs)
  x_loc <- rep(0, nc)
  x_scale <- rep(0, nc)
  j <- 0
  for (cn in covs) {
    j <- j + 1
    xx <- sub_df[[cn]]
    x_loc[j] <- mean(xx)
    x_scale[j] <- stats::sd(xx)
    x[[j]] <- (xx - x_loc[j]) / x_scale[j]
  }

  # Return
  out <- list(
    x = t(sapply(x, function(x) x)),
    x_loc = x_loc,
    x_scale = x_scale,
    nc = nc
  )
  names(out) <- paste0(names(out), "_", name)
  out
}

# Covariates
create_stan_data_covariates <- function(pd, model) {
  x_oth <- standata_scaled_covariates(pd, model, "haz")
  x_ka <- standata_scaled_covariates(pd, model, "ka")
  x_CL <- standata_scaled_covariates(pd, model, "CL")
  x_V2 <- standata_scaled_covariates(pd, model, "V2")

  # Return
  c(x_oth, x_ka, x_CL, x_V2)
}

# Create additional Stan data, interval information
create_stan_data_intervalidx <- function(t_start, t_end, t_grid, delta_grid) {
  N_int <- length(t_start)
  N_grid <- length(t_grid)
  t_start_idx_m1 <- rep(0, N_int)
  t_end_idx <- rep(0, N_int)
  cm <- rep(0, N_int)
  for (n in 1:N_int) {
    is <- 1 + sum(t_grid < t_start[n]) # start index of interval in grid
    ie <- sum(t_grid < t_end[n]) # end index of interval in grid

    # Edge case: no grid points fall within interval
    if (ie < is) {
      ie <- is
    }
    is <- min(is, N_grid)
    ie <- min(ie, N_grid)
    t_start_idx_m1[n] <- is # due to prepadded zero
    t_end_idx[n] <- ie + 1 # due to prepadded zero
    cm[n] <- (t_end[n] - t_start[n]) / ((ie - is + 1) * delta_grid)
  }

  # Return
  list(
    t_start_idx_m1 = t_start_idx_m1,
    t_end_idx = t_end_idx,
    correction_multiplier = cm
  )
}


# Average hazard per transition type, ignoring transitions that didnt
# occur
average_haz_per_ttype <- function(pd) {
  msfit <- pd$fit_mstate()
  h0 <- msfit_average_hazard(msfit) |> dplyr::arrange(.data$trans)
  df_trans <- pd$transmat$trans_df()
  df_ttype <- df_trans |> dplyr::select("trans_idx", "trans_type")
  df_ttype$trans <- df_ttype$trans_idx
  h0 <- h0 |> left_join(df_ttype, by = "trans")
  df_mean_log_h0 <- h0 |>
    dplyr::filter(avg_haz > 0) |>
    dplyr::group_by(trans_type) |>
    mutate(log_haz = log(avg_haz)) |>
    summarize(log_h0_avg = mean(log_haz)) |>
    dplyr::arrange(trans_type)
  df_ttype |> left_join(df_mean_log_h0, by = "trans_type")
}

# Edit Stan data, sort of computes which() for the binary vectors of each
# transition
which_format_for_stan <- function(x, name) {
  H <- nrow(x)
  M <- max(rowSums(x))
  N_sum <- rep(0, H)
  a <- matrix(0, H, M)
  for (h in seq_len(H)) {
    which_idx <- which(x[h, ] == 1)
    N <- length(which_idx)
    N_sum[h] <- N
    a[h, seq_len(N)] <- which_idx
  }
  out <- list(M, N_sum, a)
  names(out) <- paste0(c("D_", "sum_", "which_"), name)
  out
}


#' Create matrix indicating possible transitions given current state
#'
#' @export
#' @param legend PathData legend
#' @return a binary matrix
legend_to_PT_matrix <- function(legend) {
  trans <- legend$transition
  M <- max(c(legend$state, legend$prev_state)) # number of states
  H <- nrow(legend) # number of transitions
  checkmate::assert_true(all(trans == 1:H))
  PT <- matrix(0, H, M)
  us <- unique(legend$prev_state)
  L <- length(us)
  for (j in seq_len(L)) {
    rows <- legend |> dplyr::filter(prev_state == us[j])
    possible <- rows$transition
    PT[possible, us[j]] <- 1
  }
  colnames(PT) <- paste0("s", 1:M)
  rownames(PT) <- legend$trans_char
  PT
}


#' State/transition legend to TFI matrix
#'
#' @export
#' @param legend PathData legend
#' @return an integer matrix
legend_to_TFI_matrix <- function(legend) {
  M <- max(c(legend$state, legend$prev_state))
  TFI <- matrix(0, M, M)
  H <- nrow(legend)
  for (h in seq_len(H)) {
    TFI[legend$prev_state[h], legend$state[h]] <- legend$transition[h]
  }
  TFI
}

# Computes time since last dose for both post and pre dose measurement
time_since_last_dose <- function(t_obs_pk, last_two_times) {
  n_before_pre <- sum(last_two_times <= t_obs_pk[1])
  n_before_post <- sum(last_two_times <= t_obs_pk[2])
  if (n_before_pre == 0) {
    n_before_pre <- 1
    stop("invalid data")
  }
  if (n_before_post == 0) {
    n_before_pre <- 1
    stop("invalid data")
  }
  t_before_pre <- last_two_times[n_before_pre]
  t_before_post <- last_two_times[n_before_post]
  c(
    t_obs_pk[1] - t_before_pre,
    t_obs_pk[2] - t_before_post
  )
}
