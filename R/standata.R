# Evaluate spline basis functions and their integrals
# For each interval
# NK = number of (inner) knots
create_stan_data_spline <- function(dat, sd, k, P, NK, PT, t_max) {
  t <- dat$time
  if (is.null(t_max)) {
    t_max <- max(t)
  }
  t_pred <- seq(0, t_max, length.out = P)
  t_ev <- t[dat$transition > 0]
  knots <- place_internal_knots(t_max, NK, t_ev)
  print(knots)
  BK <- c(0, t_max)
  M <- length(knots) + k
  N <- sd$N_obs
  SBF <- array(0, dim = c(N, M))
  SBF_pred <- bspline_basis(t_pred, k, knots, BK)
  P <- length(t_pred)

  t_end <- rep(0, N)
  t_start <- rep(0, N)
  ids <- unique(dat$subject_id) # keeps order ?
  N_sub <- length(ids)
  pb <- progress::progress_bar$new(total = N_sub)
  n_obs <- 0

  # This creates also the integer subject ids (x_sub)
  x_sub <- rep(0, N)
  id_map <- NULL

  # Loop through subjects
  for (n in 1:N_sub) {
    pb$tick()
    df_n <- dat |> dplyr::filter(subject_id == ids[n])
    id_map <- rbind(id_map, data.frame(n, ids[n]))
    R <- nrow(df_n)
    for (r in 1:R) {
      n_obs <- n_obs + 1
      t_cur <- df_n$time[r]
      t_prev <- 0 # assume first interval for each subject starts from 0
      if (r > 1) {
        t_prev <- df_n$time[r - 1]
      }
      ms <- bspline_basis(t_cur, k, knots, BK)
      t_end[n_obs] <- t_cur
      t_start[n_obs] <- t_prev
      SBF[n_obs, ] <- ms
      x_sub[n_obs] <- n
    }
  }
  colnames(id_map) <- c("x_sub", "subject_id")

  # Vector indicating the start index of each subject in data
  sub_start_idx <- sapply(seq_len(N_sub), function(x) {
    # Get the index of the first occurrence of x, or NA if not found
    index <- which(x_sub == x)[1]
    if (is.na(index)) NA else index
  })

  # Grid spline evaluations
  t_grid <- seq(0.5, ceiling(t_max), by = 1) # midpoint of each grid interval
  delta_grid <- 1.0 # length of each grid interval
  N_grid <- length(t_grid)
  SBF_grid <- bspline_basis(t_grid, k, knots, BK)

  # Grid interval indices
  sd_int <- create_stan_data_intervalidx(t_start, t_end, t_grid, delta_grid)

  # Return
  sd <- list(
    N_sbf = M,
    SBF = SBF,
    SBF_pred = SBF_pred,
    SBF_grid = SBF_grid,
    N_grid = N_grid,
    t_pred = t_pred,
    t_grid = t_grid,
    delta_grid = delta_grid,
    t_end = t_end,
    t_start = t_start,
    N_pred = P,
    x_sub = x_sub,
    N_sub = max(x_sub),
    sub_start_idx = sub_start_idx
  )
  sd <- c(sd, sd_int)
  list(
    stan_data = sd,
    id_map = id_map
  )
}

# OOS id map
create_stan_data_idmap <- function(dat) {
  ids <- unique(dat$subject_id) # keeps order ?
  N_sub <- length(ids)

  # This creates also the integer subject ids (x_sub)
  x_sub <- rep(0, N_sub)
  id_map <- NULL

  # Loop through subjects
  for (n in 1:N_sub) {
    df_n <- dat |> dplyr::filter(subject_id == ids[n])
    id_map <- rbind(id_map, data.frame(n, ids[n]))
    x_sub[n] <- n
  }
  colnames(id_map) <- c("x_sub", "subject_id")

  # Vector indicating the start index of each subject in data
  sub_start_idx <- sapply(seq_len(N_sub), function(x) {
    # Get the index of the first occurrence of x, or NA if not found
    index <- which(x_sub == x)[1]
    if (is.na(index)) NA else index
  })

  # Return
  sd <- list(
    x_sub = x_sub,
    N_sub = max(x_sub),
    sub_start_idx = sub_start_idx
  )
  list(
    stan_data = sd,
    id_map = id_map
  )
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
    pk_lloq = pk_lloq
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

# Normalized covariate to stan data
standata_scaled_covariate <- function(dat, dat_oos, covariates, ssi, ssi_oos) {
  x <- list()
  x_oos <- list()
  j <- 0
  if (!is.null(ssi)) {
    dat <- dat[ssi, ]
    dat_oos <- dat[ssi_oos, ]
  }
  for (cn in covariates) {
    j <- j + 1
    xx <- dat[[cn]]
    if (is.null(xx)) {
      stop(cn, " not found in dat")
    }
    x[[j]] <- (xx - mean(xx)) / stats::sd(xx)
    x_oos[[j]] <- (dat_oos[[cn]] - mean(xx)) / stats::sd(xx)
  }

  # Return
  list(
    x = sapply(x, function(x) x),
    x_oos = sapply(x_oos, function(x) x),
    n = length(x)
  )
}

# Other covariates than sub, avg, cnt
create_stan_data_covariates <- function(dat, dat_oos, covariates,
                                        ka_covariates,
                                        CL_covariates,
                                        V2_covariates,
                                        ssi,
                                        ssi_oos) {
  covariates_rem <- setdiff(covariates, c(
    "avg", "sub", "cnt", "country_num", "country"
  ))
  x_oth <- standata_scaled_covariate(dat, dat_oos, covariates_rem, NULL, NULL)
  x_oth$x <- t(x_oth$x)
  x_oth$x_oos <- t(x_oth$x_oos)
  x_ka <- standata_scaled_covariate(dat, dat_oos, ka_covariates, ssi, ssi_oos)
  x_CL <- standata_scaled_covariate(dat, dat_oos, CL_covariates, ssi, ssi_oos)
  x_V2 <- standata_scaled_covariate(dat, dat_oos, V2_covariates, ssi, ssi_oos)
  names(x_oth) <- c("x_oth", "x_oth_oos", "N_oth")
  names(x_ka) <- c("x_ka", "x_ka_oos", "nc_ka")
  names(x_CL) <- c("x_CL", "x_CL_oos", "nc_CL")
  names(x_V2) <- c("x_V2", "x_V2_oos", "nc_V2")

  # Return
  list(
    stan_data = c(x_oth, x_ka, x_CL, x_V2),
    x_oth_names = covariates_rem
  )
}

# Create additional Stan data, interval information
create_stan_data_intervalidx <- function(t_start, t_end, t_grid, delta_grid) {
  N_obs <- length(t_start)
  N_grid <- length(t_grid)
  t_start_idx_m1 <- rep(0, N_obs)
  t_end_idx <- rep(0, N_obs)
  cm <- rep(0, N_obs)
  for (n in 1:N_obs) {
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


#' Creating Stan data list
#'
#' @export
#' @param P number of prediction points
#' @param NK number of internal spline knots
#' @param pd PD data frame
#' @param pk PK data frame
#' @param covariates hazard covariates
#' @param ka_covariates covariates that affect ka in PK model
#' @param CL_covariates covariates that affect CL in PK model
#' @param V2_covariates covariates that affect V2 in PK model
#' @param train_sub training subjects
#' @param test_sub test subjects
#' @param t_max max time of interest
#' @param do_pk should pk model be used
#' @param prior_only ignore likelihood?
#' @return A list of data for Stan.
create_stan_data <- function(pd, pk, covariates,
                             ka_covariates,
                             CL_covariates,
                             V2_covariates,
                             train_sub, test_sub,
                             P = 30, NK = 3, t_max = NULL,
                             do_pk = TRUE,
                             prior_only = FALSE) {
  checkmate::assert_logical(prior_only)
  checkmate::assert_logical(do_pk)
  checkmate::assert_class(pd, "PathData")
  checkmate::assert_integerish(P, lower = 1, len = 1)
  checkmate::assert_integerish(NK, lower = 1, len = 1)
  dt <- pd$as_transitions()

  dat <- dt$df |> dplyr::filter(subject_id %in% train_sub)
  dat_oos <- dt$df |>
    dplyr::filter(subject_id %in% test_sub) |>
    dplyr::group_by(subject_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  PT <- legend_to_PT_matrix(dt$legend)
  N_trans <- max(dt$legend$transition)
  N_obs <- nrow(dat)
  if (!is.null(pk)) {
    stopifnot(!all(duplicated(pk$subject_id)))
  }

  # Transition and at-risk indicators
  transition <- matrix(0, N_trans, N_obs)
  at_risk <- matrix(0, N_trans, N_obs)
  for (n in seq_len(N_obs)) {
    tval <- dat$transition[n]
    if (tval > 0) {
      transition[tval, n] <- 1
    }
    possible <- PT[, dat$prev_state[n]]
    at_risk[which(possible == 1), n] <- 1
  }
  ttype <- dt$legend$trans_type

  prior_cols <- which(grepl(colnames(dat), pattern = "prior_"))
  D_history <- length(prior_cols)
  history <- dat[, prior_cols]

  # Average event rates
  pd_train <- pd$filter(subject_ids_keep = train_sub)
  df_ttype <- average_haz_per_ttype(pd_train, dt)

  # Collect
  out <- list(
    mu_w0 = df_ttype$log_h0_avg,
    N_obs = N_obs,
    N_trans = N_trans,
    at_risk = at_risk,
    transition = transition,
    ttype = ttype,
    N_trans_types = max(ttype),
    N_country = max(dat$country_num),
    x_his = history,
    D_his = D_history,
    I_sub = as.numeric("sub" %in% covariates),
    I_cnt = as.numeric("cnt" %in% covariates),
    I_avg = as.numeric("avg" %in% covariates),
    x_cnt = dat$country_num,
    x_cnt_oos = dat_oos$country_num,
    x_fda = dat$first_dose_amount,
    x_fda_oos = dat_oos$first_dose_amount,
    z_mode_is = 2,
    z_mode_oos = 1
  )

  # Spline basis
  k_spline <- 3
  NK <- NK # num of internal knots
  int_dat <- create_stan_data_spline(dat, out, k_spline, P, NK, PT, t_max)
  out <- c(out, int_dat$stan_data)
  id_map_train <- int_dat$id_map

  # Oos
  idm_oos <- create_stan_data_idmap(dat_oos)
  out_add <- list(
    x_sub_oos = idm_oos$stan_data$x_sub,
    sub_start_idx_oos = idm_oos$stan_data$sub_start_idx,
    N_sub_oos = idm_oos$stan_data$N_sub
  )
  id_map_test <- idm_oos$id_map
  out <- c(out, out_add)

  # Other covariates and pk covariates
  oth <- create_stan_data_covariates(
    dat, dat_oos, covariates, ka_covariates, CL_covariates, V2_covariates,
    out$sub_start_idx,
    out$sub_start_idx_oos
  )
  out <- c(out, oth$stan_data)

  # Another way to format the transitions and at risk matrices
  # needed for vectorizing likelihood
  out <- c(out, which_format_for_stan(transition, "trans"))
  out <- c(out, which_format_for_stan(at_risk, "risk"))

  # PK data
  out <- c(out, create_stan_data_pk(pk, out, id_map_train, id_map_test))
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

# Average hazard per transition type, ignoring transitions that didnt
# occur
average_haz_per_ttype <- function(pd, dt) {
  msdata <- pd$as_msdata()$msdata
  msfit <- fit_mstate(msdata)
  h0 <- estimate_average_hazard(msfit) |> dplyr::arrange(trans)
  df_ttype <- dt$legend |> dplyr::select(transition, trans_type)
  df_ttype$trans <- df_ttype$transition
  h0 <- h0 |> left_join(df_ttype, by = "trans")
  df_mean_log_h0 <- h0 |>
    dplyr::filter(avg_haz > 0) |>
    dplyr::group_by(trans_type) |>
    mutate(log_haz = log(avg_haz)) |>
    summarize(log_h0_avg = mean(log_haz)) |>
    dplyr::arrange(trans_type)
  df_ttype <- df_ttype |> left_join(df_mean_log_h0, by = "trans_type")
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
