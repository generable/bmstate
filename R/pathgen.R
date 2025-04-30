library(progressr)
library(carrier)

# Instant hazard
evaluate_inst_hazard <- function(t, log_h0 = NULL, m_sub, t_sub, t_pred = NULL) {
  idx <- 1 + sum(t_sub <= t) # interval where t falls
  idx <- min(idx, length(m_sub)) # if t after last interval, use the last one
  if (is.null(t_pred)) {
    x <- exp(log_h0(t) + m_sub[idx])
  } else {
    x <- exp(stats::approx(t_pred, log_h0 + m_sub[idx], t)$y)
  }
  x
}

# Approximate max instant hazard on the interval [t_start, t_max]
max_inst_hazard <- function(log_h0, m_sub) {
  J <- dim(log_h0)[[2]]
  f <- rep(0, J)
  for (j in seq_len(J)) {
    f[j] <- exp(max(log_h0[1, j, ]) + max(m_sub[1, j, ]))
  }
  f
}

# Which column of a matrix is number found in
find_column_with_number <- function(matrix, number) {
  # ceiling(which(matrix == number) / nrow(matrix)) # <- slower
  for (col in 1:ncol(matrix)) {
    if (number %in% matrix[, col]) {
      return(col)
    }
  }
  stop("error in matrix")
}


# Generate transition given state and transition intensity functions
generate_transition <- function(state, log_h0 = NULL, m_sub, t_sub,
                                TFI, t_init, t_max, UB = NULL,
                                t_pred = NULL, tol = 1.03) {
  # Which transition functions are possible?
  possible <- TFI[state, ]
  possible <- possible[possible > 0]
  if (is.null(UB)) {
    UB <- max_inst_hazard(log_h0, m_sub) * tol
  }
  possible <- possible[which(UB[possible] > 1e-9)] # so rare that not even possible
  if (length(possible) == 0) {
    # Absorbing state, no transitions possible
    return(list(t = t_max, new_state = 0))
  }
  m_sub <- m_sub[1, possible, , drop = FALSE]
  UB <- UB[possible]
  if (is.null(t_pred)) {
    log_h0 <- log_h0[possible]
  } else {
    log_h0 <- log_h0[1, possible, , drop = FALSE]
  }
  # Draw event with highest upper bound first because it is most likely
  # to occur soon and be the minimum
  draw_order <- sort(UB, decreasing = TRUE, index.return = 1)$ix
  t_min_found <- t_max
  trans_idx <- 0

  # Loop through possible transitions
  for (j in draw_order) {
    if (!is.null(t_pred)) {
      h0_trans <- log_h0[1, j, ]
    } else {
      h0_trans <- log_h0[[j]]
    }
    draw <- draw_time_cnhp(
      h0_trans, m_sub[1, j, ], t_sub, UB[[j]], t_min_found, t_init, t_pred
    )
    if (draw$t < t_min_found) {
      trans_idx <- possible[j]
      t_min_found <- draw$t
    }
  }
  if (trans_idx > 0) {
    new_state <- find_column_with_number(TFI, trans_idx)
  } else {
    new_state <- 0
  }

  # Safeguard against infinite loop
  dt_min <- 1e-9
  if (t_min_found - t_init < dt_min) {
    t_min_found <- t_init + dt_min
  }

  # Return
  list(
    t = t_min_found, new_state = new_state
  )
}

# Generate path given initial state and transition intensity functions
generate_path <- function(init_state, m_sub, t_sub,
                          TFI, t_start, t_max, n_trans, UB = NULL, log_h0,
                          discretize = TRUE,
                          t_pred = NULL) {
  checkmate::assert_number(t_start, lower = 0)
  checkmate::assert_number(t_max, lower = t_start)

  # Discretize start time to start of day
  dt <- 1
  if (discretize) {
    t_start <- dt * floor(t_start / dt)
  }

  # Setup
  states <- init_state
  times <- t_start
  checkmate::assert_number(t_max, na.ok = FALSE)
  checkmate::assert_number(t_start, na.ok = FALSE)
  j <- 0

  # Generate transitions
  while (times[j + 1] < t_max) {
    j <- j + 1
    t_cur <- times[j]

    # Generate next transition
    trans <- generate_transition(
      states[j], log_h0, m_sub, t_sub, TFI, t_cur, t_max, UB, t_pred
    )
    t_next <- trans$t

    # Discretize transition time to end of day
    if (discretize) {
      t_next <- dt * ceiling(t_next / dt)
    }

    # Update time and state
    times <- c(times, t_next)
    states <- c(states, trans$new_state)
    if (j >= 100000) {
      msg <- paste0(
        "\nt = ",
        round(times[j], 2), " -> ", t_max, ": at state ", states[j], "\n ",
        "Generating a path with over ", j, " transitions, something is wrong\n"
      )
      stop(msg)
    }
  }
  L <- length(times)
  is_event <- rep(1, L)
  is_event[1] <- 0 # initial state is never an event
  is_event[L] <- 0
  states[L] <- states[L - 1]
  cbind(time = times, state = states, is_event)
}

.generate_paths_data <- function(TFI, init_state, t_start, t_max,
                                 log_h0, m_sub, t_sub,
                                 sd,
                                 sub_ids = 1, which_xsub = NULL,
                                 discretize = TRUE,
                                 n_repeats = 1,
                                 tol = 1.03, p = NULL,
                                 use_precomp = TRUE) {
  df <- NULL
  num_draws <- dim(log_h0)[1]
  stopifnot(num_draws == 1)
  if (use_precomp) {
    max_log_h0 <- apply(log_h0, 2, max, simplify = TRUE)
    log_h0 <- seq_len(sd$N_trans) |>
      purrr::map(~ stats::approxfun(sd$t_pred, log_h0[1, .x, ]))
    t_pred <- NULL
  } else {
    t_pred <- sd$t_pred
  }
  for (sub in sub_ids) {
    if (!is.null(which_xsub)) {
      subject_rows <- which_xsub[[as.character(sub)]]
      if (!is.null(t_sub)) {
        sub_t_sub <- t_sub[subject_rows]
      } else {
        sub_t_sub <- NULL
      }
      sub_m_sub <- m_sub[1, , subject_rows, drop = FALSE]
    } else {
      sub_t_sub <- t_sub
      sub_m_sub <- m_sub
    }
    if (use_precomp) {
      UB <- exp(max_log_h0 + apply(sub_m_sub, 2, max, simplify = TRUE)) * tol
    } else {
      UB <- NULL
    }
    # Generate repeat paths
    for (r in seq_len(n_repeats)) {
      path <- generate_path(
        init_state = init_state,
        log_h0 = log_h0,
        m_sub = sub_m_sub,
        t_sub = sub_t_sub,
        TFI = TFI,
        t_start = t_start,
        t_max = t_max,
        n_trans = sd$N_trans,
        UB = UB,
        discretize = discretize,
        t_pred = t_pred
      )
      df <- rbind(df, cbind(path, rep_idx = r, sub_idx = sub))
    }
  }
  if (!is.null(p)) {
    p()
  }
  as.data.frame(df)
}

# Generate multiple paths (using hazard draws after fitting a model)
generate_paths <- function(TFI, init_state, t_start, t_max,
                           log_h0, m_sub, t_sub,
                           subject_id,
                           sd, tr_names, st_names,
                           terminal_states,
                           discretize = TRUE,
                           n_repeats = 1) {
  num_draws <- dim(log_h0)[1]

  df <- seq_len(num_draws) |>
    purrr::map_dfr(
      ~ .generate_paths_data(TFI,
        init_state, t_start, t_max,
        log_h0[.x, , , drop = FALSE],
        m_sub[.x, , , drop = FALSE],
        t_sub,
        subject_id,
        sd,
        discretize = discretize,
        n_repeats = n_repeats
      ),
      .id = "draw_idx"
    ) |>
    mutate(draw_idx = as.integer(draw_idx))
  df$subject_id <- subject_id
  PathData$new(df,
    state_names = st_names,
    terminal_states = terminal_states,
    check_order = FALSE
  )
}

#' Generate all paths for some subjects
#'
#' @export
#' @param fit Stan model fit
#' @param sd Stan data
#' @param covs covariates
#' @param TFI the TFI matrix
#' @param init_state state where to initialize the paths
#' @param t_start start time
#' @param t_max_gen max time
#' @param state_names state names
#' @param trans_names transition names
#' @param terminal_states names of terminal states
#' @param df_subjects subjects data frame
#' @param num_paths number of draws for which to generate paths
#' @param oos out-of-sample mode?
#' @param discretize discretize the times?
#' @param n_repeats number of times to repeat path generation for each draw
#' @param use_future use the \code{future} package?
#' @param use_precomp use precomputation?
#' @param tol multiplier for upper bound of hazard
#' @return TODO: this function could be written more compactly
#' with less arguments
generate_paths_many_subjects <- function(fit, sd, covs,
                                         TFI, init_state, t_start, t_max_gen,
                                         state_names,
                                         trans_names,
                                         terminal_states,
                                         df_subjects,
                                         num_paths = NULL,
                                         oos = FALSE,
                                         discretize = TRUE,
                                         n_repeats = 1,
                                         use_future = TRUE,
                                         use_precomp = TRUE,
                                         tol = 1.03) {
  checkmate::assert_numeric(tol, lower = 1, len = 1)
  checkmate::assert_logical(use_future, len = 1)
  checkmate::assert_logical(use_precomp, len = 1)
  checkmate::assert_logical(discretize, len = 1)
  checkmate::assert_integerish(n_repeats, len = 1, lower = 1)
  checkmate::assert_integerish(num_paths, len = 1, lower = 1)

  # Get draws as arrays
  message(glue::glue(" * reading draws..."))
  log_h0 <- get_and_format_log_h0_draws(fit)
  m_sub <- get_and_format_log_C_haz_draws(fit, oos)

  # sampled draws based on dimension of log_h0
  num_draws <- dim(log_h0)[[1]]
  draw_ids <- seq_len(num_draws)
  if (!is.null(num_paths) && num_paths < num_draws) {
    draw_ids <- sample(draw_ids, size = num_paths, replace = FALSE)
  }
  rm(num_draws)
  sub_ids <- as.integer(df_subjects$sub_idx)

  if (isFALSE(oos)) {
    x_sub <- sd$x_sub
    t_sub <- sd$t_end
  } else {
    x_sub <- sd$x_sub_oos
    t_sub <- NULL
  }

  message(glue::glue(" * generating {length(draw_ids)*length(sub_ids)*n_repeats} paths ({length(draw_ids)*n_repeats} for each of {length(sub_ids)} subjects)..."))
  # Run path generation for the subjects (creates large data frame)

  p <- progressr::progressor(length(draw_ids))
  which_xsub <- sub_ids |>
    rlang::set_names() |>
    purrr::map(~ which(x_sub == .x))
  if (isTRUE(use_future)) {
    df <- draw_ids |>
      purrr::map_dfr(
        ~ .generate_paths_data(
          TFI = TFI, init_state = init_state, t_start = t_start, t_max = t_max_gen,
          log_h0 = log_h0[.x, , , drop = FALSE],
          m_sub = m_sub[.x, , , drop = FALSE],
          t_sub = t_sub,
          sd = sd,
          sub_ids = sub_ids,
          which_xsub = which_xsub,
          discretize = discretize,
          n_repeats = n_repeats,
          use_precomp = use_precomp,
          tol = tol, p = NULL
        ),
        .progress = TRUE,
        .id = "draw_idx",
        .options = furrr::furrr_options(
          seed = TRUE, scheduling = 1L,
          globals = c(
            "generate_transition",
            "find_column_with_number",
            "max_inst_hazard",
            "draw_time_cnhp",
            "evaluate_inst_hazard",
            ".estimate_ub",
            ".generate_paths_data",
            "generate_path"
          )
        )
      )
  } else {
    df <- draw_ids |>
      purrr::map_dfr(~ .generate_paths_data(
        TFI, init_state, t_start,
        t_max = t_max_gen,
        log_h0[.x, , , drop = FALSE],
        m_sub[.x, , , drop = FALSE],
        t_sub = t_sub,
        sd,
        sub_ids,
        which_xsub = which_xsub,
        discretize,
        n_repeats,
        use_precomp = use_precomp,
        tol, p = p
      ), .id = "draw_idx")
  }

  # Add path id & create pathdata format
  message(" * creating pathdata...")
  df$draw_idx <- as.integer(df$draw_idx)
  df$path_id <- paste0(df$sub_idx, "-", df$draw_idx, "-", df$rep_idx)
  df_link <- df |>
    dplyr::select(path_id, draw_idx, rep_idx, sub_idx) |>
    dplyr::distinct(path_id, .keep_all = TRUE) |>
    mutate(sub_idx = as.integer(sub_idx)) |>
    left_join(
      df_subjects |>
        transmute(sub_idx = as.integer(sub_idx), subject_id),
      by = join_by(sub_idx)
    )
  stopifnot(all(!is.na(df_link$subject_id)))
  df_subjects <- df_subjects |>
    semi_join(df_link, by = "subject_id")
  PathData$new(
    df_subjects, df, df_link,
    state_names = state_names,
    terminal_states = terminal_states, covs = covs,
    check_order = FALSE
  )
}


# Draw an event time from Censored Non-Homogeneous Poisson process using the
# thinning algorithm
# If drawn event time is going to be larger than t_max, then t_max is returned
draw_time_cnhp <- function(log_h0, m_sub, t_sub, lambda_ub, t_max,
                           t_init = 0, t_pred) {
  n_tries <- 0
  t <- t_init
  accepted <- FALSE
  while (!accepted) {
    n_tries <- n_tries + 1
    u1 <- stats::runif(1)
    t <- t - 1 / lambda_ub * log(u1)
    if (t > t_max) {
      return(list(t = t_max, n_tries = n_tries, censor = TRUE))
    }
    lambda_t <- evaluate_inst_hazard(t, log_h0, m_sub, t_sub, t_pred)
    if (lambda_t > lambda_ub) {
      msg <- "lambda(t) was larger than its supposed upper bound"
      msg <- paste0(
        msg, "\nt = ", t, ", t_init = ", t_init,
        ", \nlambda(t) = ", lambda_t, ", UB = ", lambda_ub
      )
      stop(msg)
    }
    p <- lambda_t / lambda_ub
    u2 <- stats::runif(1)
    if (u2 <= p) {
      accepted <- TRUE
    }
  }
  return(list(t = t, n_tries = n_tries, censor = FALSE))
}

# CF dosing
create_counterfactual_dose_data <- function(sd, dose_ss, prior_z = FALSE) {
  sd$dose_ss <- rep(dose_ss, sd$N_sub)
  sd$dose_ss_oos <- rep(dose_ss, sd$N_sub_oos)
  if (prior_z) {
    sd$z_mode_is <- 1
  } else {
    sd$z_mode_is <- 2
  }
  sd
}

predict_paths <- function(fit, stan_dat, pd, t_max_gen = NULL, oos = FALSE,
                          num_paths = NULL, sub_ids = NULL,
                          t_start = 0, init_state = 1, ...) {
  thin <- 1
  if (!is.null(num_paths)) {
    num_chains <- fit$num_chains()
    num_iter <- fit$runset$command_args() |>
      purrr::map(keep, ~ str_detect(.x, pattern = "num_samples")) |>
      purrr::map(~ str_extract(.x, "num_samples=(\\d+)", group = 1)) |>
      unlist() |>
      unique() |>
      as.integer()
    num_draws <- num_chains * num_iter
    if (num_paths < num_draws) {
      thin <- num_draws / num_paths
    }
  }
  if (!exists("mod")) {
    wdir <- tempdir(check = TRUE)
    program_file <- fs::path(wdir, "model.stan")
    brio::write_lines(fit$runset$stan_code(), program_file)
    mod <- cmdstan_model(program_file)
  }
  gq <- mod$generate_quantities(
    fit$draws() |> posterior::thin_draws(thin = thin),
    data = stan_dat$stan_data,
    parallel_chains = fit$num_chains()
  )
  TFI <- legend_to_TFI_matrix(pd$as_transitions()$legend)
  if (is.null(t_max_gen)) {
    t_max_gen <- max(stan_dat$stan_data$t_end)
  }
  if (isTRUE(oos)) {
    subject_df <- pd$subject_df |>
      inner_join(stan_dat$id_map_test, by = "subject_id") |>
      dplyr::rename(sub_idx = x_sub)
  } else {
    subject_df <- pd$subject_df |>
      inner_join(stan_dat$id_map_train, by = "subject_id") |>
      dplyr::rename(sub_idx = x_sub)
  }
  if (!is.null(sub_ids)) {
    subject_df <- subject_df |>
      dplyr::filter(subject_id %in% sub_ids)
  }

  # paths data
  paths <- generate_paths_many_subjects(
    gq, stan_dat$stan_data,
    covs = pd$covs,
    TFI, init_state, t_start, t_max_gen = t_max_gen,
    state_names = pd$state_names,
    trans_names = pd$as_transitions()$legend$trans_char,
    df_subjects = subject_df,
    terminal_states = pd$terminal_states,
    oos = oos,
    ...
  )
}

predict_cf_paths <- function(fit, stan_dat, pd, dose, t_max_gen = NULL, oos = FALSE, num_paths = NULL, sub_ids = NULL,
                             t_start = 0, init_state = 1, ...) {
  if (!is.null(num_paths)) {
    num_chains <- fit$num_chains()
    num_iter <- fit$runset$command_args() |>
      purrr::map(keep, ~ str_detect(.x, pattern = "num_samples")) |>
      purrr::map(~ str_extract(.x, "num_samples=(\\d+)", group = 1)) |>
      unlist() |>
      unique() |>
      as.integer()
    num_draws <- num_chains * num_iter
    thin <- num_draws / num_paths
  } else {
    thin <- 1
  }
  sd_cf <- create_counterfactual_dose_data(stan_dat$stan_data, dose, prior_z = FALSE)
  if (!exists("mod")) {
    wdir <- tempdir(check = TRUE)
    program_file <- fs::path(wdir, "model.stan")
    brio::write_lines(fit$runset$stan_code(), program_file)
    mod <- cmdstanr::cmdstan_model(program_file)
  }
  gq <- mod$generate_quantities(
    fit$draws() |> posterior::thin_draws(thin = thin),
    data = sd_cf,
    parallel_chains = fit$num_chains()
  )
  TFI <- legend_to_TFI_matrix(pd$as_transitions()$legend)
  if (is.null(t_max_gen)) {
    t_max_gen <- max(stan_dat$stan_data$t_end)
  }
  if (isTRUE(oos)) {
    subject_df <- pd$subject_df |>
      inner_join(stan_dat$id_map_test, by = "subject_id") |>
      mutate(cf_dose = dose) |>
      dplyr::rename(sub_idx = x_sub)
  } else {
    subject_df <- pd$subject_df |>
      inner_join(stan_dat$id_map_train, by = "subject_id") |>
      mutate(cf_dose = dose) |>
      dplyr::rename(sub_idx = x_sub)
  }
  if (!is.null(sub_ids)) {
    subject_df <- subject_df |>
      dplyr::filter(subject_id %in% sub_ids)
  }

  # paths data
  paths <- generate_paths_many_subjects(
    gq, stan_dat$stan_data,
    covs = c(pd$covs, "cf_dose"),
    TFI, init_state, t_start, t_max_gen = t_max_gen,
    state_names = pd$state_names,
    trans_names = pd$as_transitions()$legend$trans_char,
    df_subjects = subject_df,
    terminal_states = pd$terminal_states,
    oos = oos,
    ...
  )
}
