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
MultistateModelStanFit <- R6::R6Class("MultistateModelStanFit",
  private = list(
    draws = NULL,
    stan_data = NULL
  ),
  public = list(
    model = NULL,

    #' @description
    #' Create model fit object
    #'
    #' @param stan_fit A 'Stan' fit object
    #' @param stan_data The used 'Stan' data list.
    #' @param model A \code{\link{MultistateModel}}
    initialize = function(stan_fit, stan_data, model) {
      checkmate::assert_class(model, "MultistateModel")
      pars <- c(
        "weights", "log_w0", "beta_ka", "beta_V2", "beta_CL", "beta_oth",
        "sigma_pk", "log_z_pk", "log_mu_pk", "log_sig_pk", "lp__"
      )
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
  log_mu <- fit$get_draws_of("log_mu_pk")
  log_sig <- fit$get_draws_of("log_sig_pk")
  if (oos_mode) {
    log_z <- fit$get_draws_of("log_z_pk")
  } else {
    log_z <- array(0, dim = c(S, 1, sd$N_sub, sd$N_trans))
  }
  beta_ka <- fit$get_draws_of("beta_ka")
  beta_CL <- fit$get_draws_of("beta_CL")
  beta_V2 <- fit$get_draws_of("beta_V2")

  # Call exposed Stan function
  out <- list()
  for (s in seq_len(S)) {
    theta <- NULL
    if (sd$do_pk == 1) {
      theta <- compute_theta_pk(
        mat2list(log_z[s, 1, , ]),
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
#' @inheritParams msmsf_log_hazard_multipliers
msmsf_exposure <- function(fit, data = NULL) {
  # Check


  # Get draws
  beta_oth <- fit$get_draws_of("beta_oth")

  # Call exposed Stan function
  S <- fit$num_draws()
  out <- list()
  for (s in seq_len(S)) {
    out[[s]] <- compute_log_hazard_multiplier(
      sd$N_int,
      beta_oth[s, , ],
      beta_auc[s, ],
      sd$x_haz,
      x_auc[s],
      sd$ttype
    )
  }
  out
}


#' Compute log_hazard multipliers
#'
msmsf_log_hazard_multipliers <- function(fit, data = NULL) {
  # Check
  checkmate::assert_class(fit, "MultistateModelStanFit")
  if (is.null(data)) {
    sd <- fit$get_data()
  } else {
    sd <- create_stan_data(fit$model, data)
  }
  ensure_exposed_stan_functions()

  # Get draws
  beta_oth <- fit$get_draws_of("beta_oth")

  # Call exposed Stan function
  S <- fit$num_draws()
  out <- list()
  for (s in seq_len(S)) {
    out[[s]] <- compute_log_hazard_multiplier(
      sd$N_int,
      beta_oth[s, , ],
      beta_auc[s, ],
      sd$x_haz,
      x_auc[s],
      sd$ttype
    )
  }
  out
}


#' Path generation for 'MultistateModelStanFit'
#'
#' @export
#' @param fit A \code{\link{MultistateModelStanFit}} object
#' @param init_state Index of starting state
#' @param t_max Max time (start time is always 0). If \code{NULL}, the max
#' time of the model is used.
#' @param n_rep Number of repeats per draw.
#' @return A \code{\link{PathData}} object.
generate_paths <- function(fit, init_state = 1, t_max = NULL, n_rep = 10) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  checkmate::assert_integerish(n_rep, lower = 1, len = 1)

  # Get and reshape draws
  sys <- fit$model$system
  S <- fit$num_draws()
  N <- fit$stan_data$N_sub
  d <- get_inst_hazard_param_draws(fit)

  # Generate path df
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
  sub_df <- fit$data$paths$subject_df
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
#' @param fit A \code{\link{MultistateModelStanFit}} object
#' @param init_state Index of initial state
#' @param t_init Initial time
#' @param t_end End time
#' @return A data frame, where each row has the probabilities that a given subject
#' will be in each state at time \code{t_end} given that they were in
#' \code{init_state} at \code{t_init}
solve_trans_prob_fit <- function(fit, init_state = 1, t_init = 0, t_end = NULL) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  S <- fit$model$system$num_states()
  checkmate::assert_integerish(init_state, len = 1, lower = 1, upper = S)
  message("Solving transition probabilities")
  tp <- solve_trans_prob_matrix_each_subject(fit, t_init, t_end)
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
#' @param fit A \code{\link{MultistateModelStanFit}} object
#' @param t_init Initial time
#' @param t_end End time
#' @return For each subject and each draw, a matrix \code{P} where
#' \code{P[i,j]} is the probability that the system will be in state
#' \code{j} at time \code{t_end}
#' given that it is in state \code{i} at time \code{t_init}
solve_trans_prob_matrix_each_subject <- function(fit, t_init = 0, t_end = NULL) {
  checkmate::assert_class(fit, "MultistateModelStanFit")

  # Get and reshape draws
  sys <- fit$model$system
  S <- fit$num_draws()
  N <- fit$stan_data$N_sub
  d <- get_inst_hazard_param_draws(fit)
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
