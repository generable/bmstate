# Checks to be made that make sure that the model settings
# make sense for given data
prefit_checks <- function(model, data) {
  tmax <- model$get_tmax()
  pd <- data$paths
  lens <- pd$as_transitions() |>
    dplyr::mutate(time_len = .data$time - .data$time_prev) |>
    dplyr::pull(.data$time_len)
  min_len <- min(lens)
  delta_grid <- tmax / model$n_grid
  if (min_len < delta_grid) {
    msg <- paste0(
      "Shortest time interval (", min_len,
      ") is smaller than delta_grid (", delta_grid,
      "). Consider increasing n_grid or decreasing t_max of the model."
    )
    message(msg)
  }
  max_time <- max(pd$get_path_df()$time)
  if (tmax < max_time) {
    msg <- paste0(
      "Model t_max (", tmax,
      ") is smaller than max observed time in data (",
      max_time, "). Increase t_max."
    )
    stop(msg)
  }
  if (tmax > 2 * max_time) {
    warning(
      "Set model t_max is more than twice as large as largest observed time",
      " in the data. Are you sure you want to do this? If not, set smaller",
      " t_max."
    )
  }
  TRUE
}

#' Create the main 'Stan' model
#'
#' @export
#' @param ... Arguments passed to 'CmdStanR'
#' @param filepath File path to the \code{msm.stan} file.
create_stan_model <- function(filepath = NULL, ...) {
  fn <- "msm.stan"
  if (is.null(filepath)) {
    filepath <- system.file(file.path("stan", fn), package = "bmstate")
  }
  cmdstanr::cmdstan_model(filepath, ...)
}

#' Expose 'Stan' functions if they are not yet exposed
#'
#' @export
#' @param ... Passed to \code{\link{create_stan_model}}
#' @return logical telling whether they were already exposed
ensure_exposed_stan_functions <- function(...) {
  if (exists("STAN_dummy_function")) {
    return(TRUE)
  } else {
    message("Recompiling Stan model")
    mod <- create_stan_model(force_recompile = TRUE, ...)
    mod$expose_functions(global = TRUE)
  }
  FALSE
}

#' Fit a model using 'Stan'
#'
#' @description
#' \emph{NOTE:} This function has a side effect of setting normalizers.
#'
#' @export
#' @param model A \code{\link{MultistateModel}} object.
#' @param data A \code{\link{JointData}} or \code{\link{PathData}} object.
#' @param prior_only Sample from prior only?
#' @param pk_only Do not fit hazard model parameters?
#' @param return_stanfit Return also the raw 'Stan' fit object?
#' @param set_auc_normalizers Set AUC normalization based on SS doses.
#' @param filepath Passed to \code{\link{create_stan_model}}.
#' @param method Must be one of \code{"sample"} (default),
#' \code{"pathfinder"} or \code{"optimize"}.
#' @param ... Arguments passed to the \code{sample},
#' \code{pathfinder} or \code{optimize}
#' method of the 'CmdStanR' model.
#' @return A \code{\link{MultistateModelFit}} object.
fit_stan <- function(model, data, prior_only = FALSE,
                     pk_only = FALSE,
                     set_auc_normalizers = TRUE,
                     filepath = NULL,
                     return_stanfit = FALSE,
                     method = "sample", ...) {
  checkmate::assert_character(method, len = 1, min.chars = 1)
  checkmate::assert_choice(method, c("sample", "pathfinder", "optimize"))
  if (inherits(data, "PathData")) {
    # no dosing data, just use path data
    data <- JointData$new(data, NULL)
  }
  checkmate::assert_class(model, "MultistateModel")
  checkmate::assert_class(data, "JointData")
  checkmate::assert_logical(return_stanfit, len = 1)
  checkmate::assert_logical(prior_only, len = 1)
  checkmate::assert_logical(pk_only, len = 1)
  checkmate::assert_logical(set_auc_normalizers, len = 1)
  prefit_checks(model, data)

  # Get Stan model object
  stan_model <- create_stan_model(filepath = filepath)

  # Set normalizing locations and scales (side effect)
  model$set_normalizers(data)
  if (!is.null(data$dosing) && set_auc_normalizers) {
    mu_CL <- exp(-2)
    aaa <- data$dosing$dose_ss / mu_CL
    loc <- mean(aaa)
    sca <- stats::sd(aaa)
    model$set_auc_normalizers(loc, sca)
  }

  # Create Stan input list
  sd <- create_stan_data(model, data, prior_only, pk_only)

  # Call 'Stan'
  if (method == "pathfinder") {
    stan_fit <- stan_model$pathfinder(data = sd, ...)
  } else if (method == "optimize") {
    stan_fit <- stan_model$optimize(data = sd, ...)
  } else {
    stan_fit <- stan_model$sample(data = sd, ...)
  }

  # Return
  pars <- c(
    "weights", "log_w0", "beta_ka", "beta_V2", "beta_CL", "beta_oth",
    "beta_auc", "sigma_pk", "log_z_pk", "log_mu_pk", "log_sig_pk", "lp__"
  )
  draws <- create_rv_list(stan_fit, pars)
  diag <- NULL
  if (method == "sample") {
    diag <- stan_fit$diagnostic_summary()
  }
  sfun <- function(x) x$summary()
  quiet_sumr <- purrr::quietly(sfun)
  smr <- quiet_sumr(stan_fit)
  info <- list(
    diag = diag,
    runset = stan_fit$runset,
    summary = smr$result
  )
  fit <- MultistateModelFit$new(data, sd, model, draws, info)
  if (!return_stanfit) {
    return(fit)
  }
  list(
    fit = fit,
    stanfit = stan_fit
  )
}


#' Creating Stan data list
#'
#' @export
#' @inheritParams fit_stan
#' @return A list of data for Stan.
create_stan_data <- function(model, data, prior_only = FALSE, pk_only = FALSE) {
  checkmate::assert_class(model, "MultistateModel")
  checkmate::assert_class(data, "JointData")
  pd <- data$paths
  tm <- pd$transmat
  check_equal_transmats(tm, model$system$tm())
  checkmate::assert_logical(prior_only, len = 1)
  delta_grid <- model$system$get_tmax() / model$n_grid

  # Initial Stan data
  stan_dat <- c(
    create_stan_data_idx_sub(pd),
    create_stan_data_transitions(pd),
    create_stan_data_spline(pd, model, delta_grid),
    create_stan_data_covariates(pd, model),
    create_stan_data_trans_types(pd),
    create_stan_data_pk(data, model)
  )

  # Likelihood flags
  flags <- list(
    omit_lik_haz = as.integer(prior_only),
    omit_lik_pk = as.integer(prior_only),
    do_pk = as.integer(model$has_pk()),
    do_haz = as.integer(!pk_only)
  )

  # Return
  c(stan_dat, flags)
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
create_stan_data_transitions <- function(pd) {
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
    at_risk[tm$possible_transitions_from(dat$from[n]), n] <- 1
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

# Map of subject id to integer
stan_data_id_map <- function(pd) {
  checkmate::assert_class(pd, "PathData")
  sid <- pd$unique_subjects()
  idx <- seq_len(length(sid))
  data.frame(subject_index = idx, subject_id = sid)
}

# Which subject does each interval correspond to
create_stan_data_idx_sub <- function(pd) {
  im <- stan_data_id_map(pd)
  dt <- pd$as_transitions() |> dplyr::left_join(im, by = "subject_id")
  idx_sub <- as.integer(dt$subject_index)
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
  if (delta_grid > 0.25 * t_max) {
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

# PK data
create_stan_data_pk <- function(data, model) {
  checkmate::assert_class(data, "JointData")
  if (model$has_pk()) {
    pk_obs <- data$paths$subject_df |>
      dplyr::select("t_pre", "t_post", "conc_pre", "conc_post", "pk_lloq")
    pk_obs <- as.matrix(pk_obs)
    last_two_times <- t(sapply(data$dosing$times, function(x) x))
    last_two_doses <- t(sapply(data$dosing$doses, function(x) x))
    dose_ss <- data$dosing$dose_ss
    tau_ss <- data$dosing$tau_ss
  } else {
    N_sub <- nrow(data$paths$subject_df)
    pk_obs <- matrix(1, N_sub, 5)
    last_two_doses <- matrix(0, N_sub, 2)
    last_two_times <- matrix(1, N_sub, 2)
    dose_ss <- rep(0, N_sub)
    tau_ss <- 0
  }
  t_obs_pk <- pk_obs[, 1:2]
  conc_pk <- pk_obs[, 3:4]
  pk_lloq <- as.numeric(pk_obs[, 5])
  an <- model$get_auc_normalizers()

  # Return
  sd <- list(
    t_obs_pk = t_obs_pk,
    conc_pk = conc_pk,
    dose_ss = dose_ss,
    last_two_times = last_two_times,
    last_two_doses = last_two_doses,
    pk_lloq = pk_lloq,
    I_auc = as.numeric(model$has_pk()),
    tau_ss = tau_ss,
    auc_loc = an$loc,
    auc_scale = an$scale
  )
  c(sd, create_stan_data_time_since_last_pk(sd))
}

# Time since last dose
create_stan_data_time_since_last_pk <- function(sd) {
  N_sub <- length(sd$dose_ss)
  t_since_last_pk <- matrix(0, N_sub, 2)
  for (n in seq_len(N_sub)) {
    tsl <- time_since_last_dose(
      as.numeric(sd$t_obs_pk[n, ]),
      as.numeric(sd$last_two_times[n, ])
    )
    t_since_last_pk[n, ] <- tsl
  }
  list(
    t_since_last_pk = t_since_last_pk
  )
}

# Normalized covariates to Stan data
standata_scaled_covariates <- function(pd, model, name) {
  covs <- model$data_covs(name)
  sub_df <- pd$subject_df
  N <- nrow(sub_df)
  x <- list()
  nc <- length(covs)
  j <- 0
  norms <- model$get_normalizers()
  for (cn in covs) {
    j <- j + 1
    xx <- sub_df[[cn]]
    xj_loc <- norms$locations[[cn]]
    xj_scale <- norms$scales[[cn]]
    if (is.null(xj_loc) || is.null(xj_scale)) {
      stop("no normalizers set for ", cn)
    }
    x[[j]] <- (xx - xj_loc) / xj_scale
  }
  x <- sapply(x, function(x) x)
  if (length(dim(x)) != 2) {
    x <- array(dim = c(N, 0))
  }
  if (nrow(x) != N) {
    stop("error in standata_scaled_covariates")
  }

  # Return
  out <- list(
    x = x,
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
  out <- c(x_oth, x_ka, x_CL, x_V2)
  out$x_haz <- t(out$x_haz)
  out
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

# Average hazard per transition type, ignoring transitions that did not
# occur
average_haz_per_ttype <- function(pd) {
  msfit <- pd$fit_mstate()
  h0 <- msfit_average_hazard(msfit) |> dplyr::arrange(.data$trans)
  df_trans <- pd$transmat$trans_df()
  df_ttype <- df_trans |> dplyr::select("trans_idx", "trans_type")
  df_ttype$trans <- df_ttype$trans_idx
  h0 <- h0 |> left_join(df_ttype, by = "trans")
  df_mean_log_h0 <- h0 |>
    dplyr::filter(.data$avg_haz > 0) |>
    dplyr::group_by(.data$trans_type) |>
    mutate(log_haz = log(.data$avg_haz)) |>
    summarize(log_h0_avg = mean(.data$log_haz)) |>
    dplyr::arrange(.data$trans_type)
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
