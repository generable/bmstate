#' Solve the Kolmogorov forward equation of a Markovian multistate system
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param t A time grid
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_w0 A vector of length \code{n_trans}
#' @param log_m A vector of length \code{n_trans}
#' @return Value given by \code{deSolve::ode}.
solve_time_evolution <- function(system, t, log_w0, w = NULL, log_m = NULL) {
  checkmate::assert_class(system, "MultistateSystem")
  S <- system$num_states()
  H <- system$num_trans()
  W <- system$num_weights()
  if (is.null(w)) {
    w <- matrix(0, H, W)
  }
  if (is.null(log_m)) {
    log_m <- rep(0, H)
  }
  checkmate::assert_numeric(t, min.len = 1)
  checkmate::assert_matrix(w, ncols = W, nrows = H)
  checkmate::assert_numeric(log_w0, len = H)
  checkmate::assert_numeric(log_m, len = H)

  odefun <- function(time, y, parms) {
    P <- matrix(y, S, S)
    Lambda <- system$intensity_matrix(time, log_w0, w, log_m)
    dydt <- as.vector(P %*% Lambda)
    list(dydt)
  }
  P0 <- diag(1, S, S)
  y0 <- as.vector(P0)
  deSolve::ode(y0, t, odefun, NULL, method = "ode45")
}

#' Solve the transition probability matrix
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param log_w0 A vector of length \code{n_trans}
#' @param t_start Initial time
#' @param t_end End time
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_m A vector of length \code{n_trans}
#' @return A matrix \code{P} where \code{P[i,j]} is the probability that
#' the system will be in state \code{j} at time \code{t_end} given that it
#' is in state \code{i} at time \code{t_start}
solve_trans_prob_matrix <- function(system, log_w0, w = NULL,
                                    log_m = NULL, t_start = 0, t_end = NULL) {
  checkmate::assert_class(system, "MultistateSystem")
  if (is.null(t_end)) {
    t_end <- system$get_tmax()
  }
  checkmate::assert_number(t_start, lower = 0)
  checkmate::assert_number(t_end, lower = t_start + 1e-9)
  H <- system$num_trans()
  W <- system$num_weights()
  S <- system$num_states()
  if (is.null(w)) {
    w <- matrix(0, H, W)
  }
  if (is.null(log_m)) {
    log_m <- rep(0, H)
  }
  t <- c(t_start, t_end)
  kfe <- solve_time_evolution(system, t, log_w0, w, log_m)
  P <- matrix(kfe[2, 2:ncol(kfe)], S, S)
  cn <- system$tm()$states
  colnames(P) <- cn
  rownames(P) <- cn
  P
}


#' Solve transition probabilities for each subject using 'MultistateModelFit'
#'
#' @export
#' @inheritParams msmfit_pk_params
#' @param t_start Initial time.
#' @param t_end End time.
#' @return A data frame, where each row has the probabilities that a given subject
#' will be in each state at time \code{t_end} given that they were in
#' \code{init_state} at \code{t_start}
solve_trans_prob_fit <- function(fit, t_start = 0, t_end = NULL,
                                 data = NULL) {
  checkmate::assert_class(fit, "MultistateModelFit")
  message("Solving transition probabilities")
  tp <- solve_trans_prob_matrix_each_subject(fit, t_start, t_end, data = data)
  message("Formatting")

  init_states <- msmfit_state_at(t_start, fit, data) |>
    dplyr::select("state", "subject_id")
  N_sub <- nrow(init_states)
  S <- fit$model$system$num_states()

  df <- matrix(0, N_sub, 1 + S)
  for (j in seq_len(N_sub)) {
    init_state_j <- init_states$state[j]
    sid <- init_states$subject_id[j]
    inds <- which(tp$index_df$subject_id == sid)
    Pi <- tp$P[inds, , ]
    if (length(dim(Pi)) > 2) {
      Pi <- apply(Pi, c(2, 3), mean)
    }
    Pi <- Pi[init_state_j, ]
    df[j, ] <- c(sid, Pi)
  }
  df <- data.frame(df)
  colnames(df) <- c("subject_id", fit$model$system$tm()$states)
  df
}

#' Solve transition probabilities for each subject in 'MultistateModelFit'
#'
#' @inheritParams solve_trans_prob_fit
#' @return For each subject and each draw, a matrix \code{P} where
#' \code{P[i,j]} is the probability that the system will be in state
#' \code{j} at time \code{t_end}
#' given that it is in state \code{i} at time \code{t_start}
solve_trans_prob_matrix_each_subject <- function(fit, t_start = 0, t_end = NULL,
                                                 data = NULL) {
  checkmate::assert_class(fit, "MultistateModelFit")
  checkmate::assert_number(t_start, lower = 0)

  # Get and reshape draws
  sys <- fit$model$system
  S <- fit$num_draws()
  sd <- msmfit_stan_data(fit, data)
  N <- sd$N_sub
  d <- msmfit_inst_hazard_param_draws(fit, data)
  NS <- nrow(d$df)
  pb <- progress::progress_bar$new(total = NS)
  K <- sys$num_states()
  A <- array(0, dim = c(NS, K, K))
  for (j in seq_len(NS)) {
    pb$tick()
    wj <- d$w[j, , ]
    if (is.null(dim(wj))) {
      wj <- matrix(wj, 1, length(wj))
    }
    A[j, , ] <- solve_trans_prob_matrix(
      sys, d$log_w0[j, ], wj, d$log_m[j, ], t_start, t_end
    )
  }
  list(P = A, index_df = d$df)
}
