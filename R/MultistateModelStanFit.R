# Helper
create_rv_list <- function(stan_fit, names) {
  out <- list()
  j <- 0
  names_out <- NULL
  for (name in names) {
    tryCatch(
      {
        rv <- rv(stan_fit, name)
        j <- j + 1
        out[[j]] <- rv
        names_out <- c(names_out, name)
      },
      error = function(e) {
      }
    )
  }
  names(out) <- names_out
  out
}

#' Minimal fit class
#'
#' @export
#' @field model The \code{\link{MultistateModel}}
#' @field data A  \code{\link{JointData}} object
MultistateModelStanFit <- R6::R6Class("MultistateModelStanFit",
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
    #' @param stan_fit A 'Stan' fit object
    #' @param stan_data The used 'Stan' data list.
    #' @param model A \code{\link{MultistateModel}}
    initialize = function(data, stan_fit, stan_data, model) {
      checkmate::assert_class(data, "JointData")
      checkmate::assert_class(model, "MultistateModel")
      pars <- c(
        "weights", "log_w0", "beta_ka", "beta_V2", "beta_CL", "beta_oth",
        "beta_auc",
        "sigma_pk", "log_z_pk", "log_mu_pk", "log_sig_pk", "lp__"
      )
      self$data <- data
      self$model <- model
      private$draws <- create_rv_list(stan_fit, pars)
      private$stan_data <- stan_data
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
      S <- self$num_draws()
      x1 <- paste0("A MultistateModelStanFit with ", S, " draws.")
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

    #' Plot baseline hazard distribution
    #'
    #' @param t times where to evaluate the baseline hazards
    #' @param ci_alpha width of credible interval
    plot_h0 = function(t = NULL, ci_alpha = 0.95) {
      df <- self$h0_dist(t, ci_alpha)
      legend <- self$model$system$tm()$trans_df()
      df <- df |> dplyr::left_join(legend, by = "trans_idx")
      df$trans_type <- as.factor(df$trans_type)
      capt <- paste0("Median and ", 100 * ci_alpha, "% CrI")
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
  checkmate::assert_class(fit, "MultistateModelStanFit")
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
#' @param fit A \code{\link{MultistateModelStanFit}} object
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
  if (sd$do_pk == 1) {
    beta_auc <- fit$get_draws_of("beta_auc")
  } else {
    beta_auc <- array(0, dim = c(S, 0, sd$N_trans_types))
  }

  # Call exposed Stan function for each draw (not optimal)
  out <- list()
  N_sub <- sd$N_sub
  first_indices <- sapply(seq_len(N_sub), function(x) which(sd$idx_sub == x)[1])
  for (s in seq_len(S)) {
    if (sd$do_pk == 1) {
      ba <- list(beta_auc[s, 1, ])
      aa <- list(auc[[s]][sd$idx_sub])
    } else {
      ba <- NULL
      aa <- NULL
    }
    out[[s]] <- compute_log_hazard_multiplier(
      sd$N_int,
      mat2list(t(beta_oth[s, , ])),
      ba,
      mat2list(t(sd$x_haz[, sd$idx_sub])),
      aa,
      sd$ttype
    )
    out[[s]] <- out[[s]][first_indices, ]
  }
  out
}


# Log baseline hazard distribution at times t
msmsf_log_baseline_hazard <- function(fit, t = NULL) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
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


#' Path generation for 'MultistateModelStanFit'
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
  checkmate::assert_class(fit, "MultistateModelStanFit")
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

#' Solve transition probabilities for each subject in 'MultistateModelStanFit'
#'
#' @export
#' @inheritParams msmsf_pk_params
#' @param init_state Index of initial state
#' @param t_init Initial time
#' @param t_end End time
#' @return A data frame, where each row has the probabilities that a given subject
#' will be in each state at time \code{t_end} given that they were in
#' \code{init_state} at \code{t_init}
solve_trans_prob_fit <- function(fit, init_state = 1, t_init = 0, t_end = NULL,
                                 data = NULL) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  S <- fit$model$system$num_states()
  checkmate::assert_integerish(init_state, len = 1, lower = 1, upper = S)
  message("Solving transition probabilities")
  tp <- solve_trans_prob_matrix_each_subject(fit, t_init, t_end, data = data)
  message("Formatting")
  NS <- length(tp$subject_index)
  us <- unique(tp$subject_index)
  N_sub <- max(us)
  df <- matrix(0, N_sub, 1 + S)
  j <- 0
  for (sidx in us) {
    j <- j + 1
    inds <- which(tp$subject_index == sidx)
    Pi <- apply(tp$P[inds, , ], c(2, 3), mean)[init_state, ]
    df[j, ] <- c(sidx, Pi)
  }
  df <- data.frame(df)
  colnames(df) <- c("subject_index", fit$model$system$tm()$states)
  df
}

#' Solve transition probabilities for each subject in 'MultistateModelStanFit'
#'
#' @inheritParams msmsf_pk_params
#' @param t_init Initial time
#' @param t_end End time
#' @return For each subject and each draw, a matrix \code{P} where
#' \code{P[i,j]} is the probability that the system will be in state
#' \code{j} at time \code{t_end}
#' given that it is in state \code{i} at time \code{t_init}
solve_trans_prob_matrix_each_subject <- function(fit, t_init = 0, t_end = NULL,
                                                 data = NULL) {
  checkmate::assert_class(fit, "MultistateModelStanFit")

  # Get and reshape draws
  sys <- fit$model$system
  S <- fit$num_draws()
  sd <- msmsf_stan_data(fit, data)
  N <- sd$N_sub
  d <- msmsf_inst_hazard_param_draws(fit, data)
  NS <- length(d$subject_index)
  pb <- progress::progress_bar$new(total = NS)
  K <- sys$num_states()
  A <- array(0, dim = c(NS, K, K))
  for (j in seq_len(NS)) {
    pb$tick()
    A[j, , ] <- solve_trans_prob_matrix(
      sys, d$log_w0[j, ], d$w[j, , ], d$log_m[j, ], t_init, t_end
    )
  }
  c(list(P = A), d[c("subject_index", "draw_index")])
}
