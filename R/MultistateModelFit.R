#' Minimal fit class
#'
#' @export
#' @field model The \code{\link{MultistateModel}}
#' @field data A  \code{\link{JointData}} object
MultistateModelFit <- R6::R6Class("MultistateModelFit",
  private = list(
    draws = NULL,
    stan_data = NULL
  ),
  public = list(
    model = NULL,
    data = NULL,

    #' @description
    #' Create model fit object
    #'
    #' @param data Data used to create the fit.
    #' @param draws A named list of rvars.
    #' @param stan_data The used 'Stan' data list.
    #' @param model A \code{\link{MultistateModel}}
    initialize = function(data, stan_data, model, draws) {
      checkmate::assert_class(data, "JointData")
      checkmate::assert_class(model, "MultistateModel")
      self$data <- data
      self$model <- model
      private$draws <- draws
      private$stan_data <- stan_data
    },

    #' @description Create a reduced version with only a single draw, corresponding
    #' to the mean of original draws.
    #' @return A new \code{\link{MultistateModelFit}} object.
    mean_fit = function() {
      draws <- lapply(private$draws, rvar_to_mean_rvar)
      MultistateModelFit$new(self$data, private$stan_data, self$model, draws)
    },

    #' @description Check if fit is a point estimate
    #' @return a logical value
    is_point_estimate = function() {
      self$num_draws() == 1
    },

    #' @description Extract data list
    #'
    get_data = function() {
      private$stan_data
    },

    #' @description Extract draws as \code{rvar}s
    #'
    #' @param name Param/quantity name
    get_draws = function(name = NULL) {
      private$draws[[name]]
    },

    #' Draws in a raw array with same shape as Stan variable
    #'
    #' @param name Param/quantity name of \code{x}
    #' @return Array with dimension \code{c(ndraws(x), dim(x))}
    get_draws_of = function(name) {
      posterior::draws_of(self$get_draws(name), with_chains = FALSE)
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      if (self$is_point_estimate()) {
        x1 <- paste0("A point-estimate MultistateModelFit.")
      } else {
        S <- self$num_draws()
        x1 <- paste0("A MultistateModelFit with ", S, " draws.")
      }
      msg <- paste(x1, "\n", sep = "\n")
      cat(msg)
    },

    #' @description Plot used basis functions (grid)
    plot_basisfun = function() {
      sd <- self$get_data()
      t <- rep(sd$t_grid, sd$N_sbf)
      y <- as.vector(sd$SBF_grid)
      idx <- as.factor(rep(1:sd$N_sbf, each = sd$N_grid))
      df <- data.frame(t, y, idx)
      ggplot(df, aes(.data$t, .data$y, color = .data$idx)) +
        geom_line() +
        ggtitle("Basis functions")
    },

    #' @description Plot PK fit.
    #' @param max_num_subjects Max number of subjects to show.
    #' @param data Data for which to predict the concentration. If \code{NULL},
    #' training data is used.
    #' @param L number of grid points for each subject
    #' @param timescale scale of time
    #' @param n_prev number of previous doses to show fit for
    plot_pk = function(max_num_subjects = 12, data = NULL, L = 30,
                       timescale = 24, n_prev = 3) {
      checkmate::assert_integerish(L, len = 1)
      if (is.null(data)) {
        data <- self$data
      }
      pkpar <- msmsf_pk_params(self, data = data)
      theta <- pkpar[[1]]
      trange <- sapply(data$dosing$times, range)
      N <- nrow(theta)
      ts <- list()
      for (j in 1:N) {
        ts[[j]] <- seq(
          trange[1, j] - timescale * n_prev,
          trange[2, j],
          length.out = L
        ) + timescale
      }
      pksim <- data$dosing$simulate_pk(ts, theta)
      pltd <- data$plot_dosing(df_fit = pksim, max_num_subjects = max_num_subjects)
      pltd
    },

    #' Plot baseline hazard distribution
    #'
    #' @param t times where to evaluate the baseline hazards
    #' @param ci_alpha width of credible interval
    plot_h0 = function(t = NULL, ci_alpha = 0.95) {
      df <- self$h0_dist(t, ci_alpha)
      legend <- self$model$system$tm()$trans_df()
      df <- df |> dplyr::left_join(legend, by = "trans_idx")
      df$trans_type <- as.factor(df$trans_type)
      if (!self$is_point_estimate()) {
        capt <- paste0("Median and ", 100 * ci_alpha, "% CrI")
      } else {
        capt <- paste0("Point-estimate fit")
      }

      ggplot(df, aes(
        x = .data$time, y = .data$median, ymin = .data$lower,
        ymax = .data$upper, color = .data$trans_type, fill = .data$trans_type
      )) +
        geom_ribbon(alpha = 0.7) +
        geom_line() +
        facet_wrap(. ~ .data$trans_char) +
        ylab("Baseline hazard") +
        theme(legend.position = "none") +
        labs(caption = capt) +
        scale_y_log10()
    },

    #' Baseline hazard distribution
    #'
    #' @param t times where to evaluate the baseline hazards
    #' @param ci_alpha width of credible interval
    h0_dist = function(t = NULL, ci_alpha = 0.95) {
      checkmate::assert_number(ci_alpha, lower = 0, upper = 1)
      LB <- (1 - ci_alpha) / 2
      UB <- 1 - LB
      df <- msmsf_log_baseline_hazard(self, t)
      df |>
        dplyr::group_by(.data$time, .data$trans_idx) |>
        dplyr::summarize(
          median = stats::median(exp(.data$log_h0)),
          upper = stats::quantile(exp(.data$log_h0), UB),
          lower = stats::quantile(exp(.data$log_h0), LB)
        ) |>
        dplyr::ungroup()
    },

    #' @description
    #' Full names of parameters that start with \code{log_z_}.
    log_z_pars = function() {
      nams <- names(self$draws())
      match <- grepl(nams, pattern = "log_z_")
      nams[which(match)]
    },

    #' @description Get number of draws
    num_draws = function() {
      posterior::ndraws(self$get_draws("lp__"))
    }
  )
)

# Helper
msmsf_stan_data <- function(fit, data = NULL) {
  checkmate::assert_class(fit, "MultistateModelFit")
  if (is.null(data)) {
    sd <- fit$get_data()
  } else {
    sd <- create_stan_data(fit$model, data)
  }
  ensure_exposed_stan_functions()
  sd
}

# Helper
mat2list <- function(mat) {
  lapply(seq_len(ncol(mat)), function(j) mat[, j])
}


#' Evaluate PK parameters
#' @export
#' @param fit A \code{\link{MultistateModelFit}} object
#' @param data A \code{\link{JointData}} object. If \code{NULL}, the
#' data used to fit the model is used. If not \code{NULL}, out-of-sample
#' mode is assumed (new subjects).
#' @return A list with length equal to number of draws.
msmsf_pk_params <- function(fit, data = NULL) {
  oos_mode <- is.null(data)
  sd <- msmsf_stan_data(fit, data)
  S <- fit$num_draws()

  # Extract
  log_mu <- fit$get_draws_of("log_mu_pk")
  log_sig <- fit$get_draws_of("log_sig_pk")
  if (oos_mode) {
    log_z <- fit$get_draws_of("log_z_pk")
  } else {
    log_z <- array(0, dim = c(S, 1, sd$N_sub, 3))
  }
  beta_ka <- fit$get_draws_of("beta_ka")
  beta_CL <- fit$get_draws_of("beta_CL")
  beta_V2 <- fit$get_draws_of("beta_V2")

  # Call exposed Stan function for each draw (not optimal)
  out <- list()
  for (s in seq_len(S)) {
    theta <- NULL
    if (sd$do_pk == 1) {
      theta <- compute_theta_pk(
        mat2list(t(log_z[s, 1, , ])),
        log_mu[s, 1, ],
        log_sig[s, 1, ],
        beta_ka[s, 1, ],
        beta_CL[s, 1, ],
        beta_V2[s, 1, ],
        mat2list(t(sd$x_ka)),
        mat2list(t(sd$x_CL)),
        mat2list(t(sd$x_V2))
      )
    }
    out[[s]] <- theta
  }
  out
}

#' Compute exposure
#'
#' @export
#' @inheritParams msmsf_pk_params
msmsf_exposure <- function(fit, data = NULL) {
  # Get draws
  sd <- msmsf_stan_data(fit, data)
  pkpar <- msmsf_pk_params(fit, data)

  # Call exposed Stan function
  S <- fit$num_draws()
  out <- list()
  for (s in seq_len(S)) {
    if (sd$do_pk == 1) {
      x_auc <- sd$dose_ss / pkpar[[s]][, 2] # D/CL
    } else {
      x_auc <- NULL
    }
    out[[s]] <- x_auc
  }
  out
}


#' Compute log_hazard multipliers
#'
#' @export
#' @inheritParams msmsf_pk_params
#' @return A list of length \code{n_draws} where each element is a
#' matrix of shape \code{n_subject} x \code{n_transitions}
msmsf_log_hazard_multipliers <- function(fit, data = NULL) {
  # Get draws
  sd <- msmsf_stan_data(fit, data)
  auc <- msmsf_exposure(fit, data)
  S <- fit$num_draws()
  beta_oth <- fit$get_draws_of("beta_oth")
  if (is.null(beta_oth)) {
    beta_oth <- array(0, dim = c(S, 0, sd$N_trans_types))
  }
  if (sd$do_pk == 1) {
    beta_auc <- fit$get_draws_of("beta_auc")
  } else {
    beta_auc <- array(0, dim = c(S, 0, sd$N_trans_types))
  }

  # Call exposed Stan function for each draw (not optimal)
  out <- list()
  N_sub <- sd$N_sub
  first_indices <- sapply(seq_len(N_sub), function(x) which(sd$idx_sub == x)[1])
  if (sd$nc_haz > 0) {
    x_haz_long <- sd$x_haz[, sd$idx_sub]
  } else {
    x_haz_long <- array(0, dim = c(0, sd$N_int))
  }

  for (s in seq_len(S)) {
    if (sd$do_pk == 1) {
      ba <- list(beta_auc[s, 1, ])
      aa <- list(auc[[s]][sd$idx_sub])
    } else {
      ba <- NULL
      aa <- NULL
    }
    if (sd$nc_haz == 0) {
      r <- matrix(0, sd$N_sub, sd$N_trans)
    } else {
      r <- compute_log_hazard_multiplier(
        sd$N_int,
        mat2list(t(beta_oth[s, , ])),
        ba,
        mat2list(t(x_haz_long)),
        aa,
        sd$ttype
      )
      r <- r[first_indices, , drop = FALSE]
    }
    out[[s]] <- r
  }
  out
}


# Log baseline hazard distribution at times t
msmsf_log_baseline_hazard <- function(fit, t = NULL) {
  checkmate::assert_class(fit, "MultistateModelFit")
  sys <- fit$model$system
  if (is.null(t)) {
    t <- seq(0, sys$get_tmax(), length.out = 30)
  }
  checkmate::assert_numeric(t, min.len = 2)
  SBF <- sys$basisfun_matrix(t)
  w <- fit$get_draws_of("weights") # dim = c(S, H, W)
  log_w0 <- fit$get_draws_of("log_w0") # dim = c(S, H)
  S <- fit$num_draws()
  N <- length(t)
  H <- sys$num_trans()
  log_h0 <- NULL
  trans_idx <- NULL
  draw_idx <- NULL
  time <- NULL
  for (h in seq_len(H)) {
    for (s in seq_len(S)) {
      log_h0_s <- sys$log_baseline_hazard(NULL, log_w0[s, h], w[s, h, ], SBF)
      log_h0 <- c(log_h0, log_h0_s)
      time <- c(time, t)
      draw_idx <- c(draw_idx, rep(s, N))
      trans_idx <- c(trans_idx, rep(h, N))
    }
  }
  data.frame(draw_idx, trans_idx, log_h0, time)
}


#' Extract and reshape draws of instant hazard related parameters
#'
#' @export
#' @inheritParams msmsf_pk_params
#' @return a list with elements \code{log_m}, \code{log_w0}, \code{w}, each
#' of which is an array where the first dimension is number of subjects times
#' number of draws
msmsf_inst_hazard_param_draws <- function(fit, data = NULL) {
  sd <- msmsf_stan_data(fit, data)
  log_m <- msmsf_log_hazard_multipliers(fit, data)
  S <- fit$num_draws()
  N <- sd$N_sub
  w <- fit$get_draws_of("weights")
  log_w0 <- fit$get_draws_of("log_w0")
  w_rep <- abind::abind(replicate(N, w, simplify = FALSE), along = 1)
  log_w0_rep <- abind::abind(replicate(N, log_w0, simplify = FALSE), along = 1)
  log_m_reshaped <- do.call(rbind, log_m)
  list(
    log_m = log_m_reshaped,
    log_w0 = log_w0_rep,
    w = w_rep,
    subject_index = rep(seq_len(N), times = S),
    draw_index = rep(seq_len(S), each = N)
  )
}


#' Path generation for 'MultistateModelFit'
#'
#' @export
#' @inheritParams msmsf_pk_params
#' @param init_state Index of starting state
#' @param t_max Max time (start time is always 0). If \code{NULL}, the max
#' time of the model is used.
#' @param n_rep Number of repeats per draw.
#' @return A \code{\link{PathData}} object.
generate_paths <- function(fit, init_state = 1, t_max = NULL, n_rep = 10,
                           data = NULL) {
  checkmate::assert_class(fit, "MultistateModelFit")
  checkmate::assert_integerish(n_rep, lower = 1, len = 1)
  sd <- msmsf_stan_data(fit, data)
  log_m <- msmsf_log_hazard_multipliers(fit, data)

  # Get and reshape draws
  sys <- fit$model$system
  S <- fit$num_draws()
  N <- sd$N_sub
  message("Computing hazard multipliers")
  d <- msmsf_inst_hazard_param_draws(fit, data)

  # Generate path df
  message("Generating paths")
  path_df <- sys$simulate(d$w, d$log_w0, d$log_m, init_state, t_max, n_rep)

  # Create indices for link
  subject_index <- rep(d$subject_index, times = n_rep)
  draw_index <- rep(d$draw_index, times = n_rep)
  rep_index <- rep(seq_len(n_rep), each = N * S)

  # Check
  stopifnot(length(subject_index) == N * S * n_rep)
  stopifnot(length(draw_index) == N * S * n_rep)
  stopifnot(length(rep_index) == N * S * n_rep)

  # Create link df
  if (is.null(data)) {
    dat <- fit$data
  } else {
    dat <- data
  }
  sub_df <- dat$paths$subject_df
  link_df <- data.frame(
    path_id = seq_len(N * S * n_rep),
    subject_id = sub_df$subject_id[subject_index],
    draw_idx = draw_index,
    rep_idx = rep_index
  )

  # Create PathData object
  PathData$new(
    subject_df = sub_df,
    path_df = path_df,
    link_df = link_df,
    transmat = fit$model$system$tm(),
    covs = fit$model$data_covs()
  )
}
