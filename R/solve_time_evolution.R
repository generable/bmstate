#' Solve the Kolmogorov forward equation of a Markovian multistate system
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param t A time grid
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_w0 A vector of length \code{n_trans}
#' @param log_m A vector of length \code{n_trans}
#' @param ... Arguments passed to \code{deSolve::ode()}.
#' @return Value given by \code{deSolve::ode()}.
solve_time_evolution <- function(system, t, log_w0, w = NULL, log_m = NULL,
                                 ...) {
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
  deSolve::ode(y0, t, odefun, NULL, method = "ode45", ...)
}

#' Solve the transition probability matrix for each given output time
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param log_w0 A vector of length \code{n_trans}
#' @param t_start Initial time
#' @param t_out End times (vector)
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_m A vector of length \code{n_trans}
#' @param ... Arguments passed to \code{deSolve::ode()}.
#' @return An array \code{P} where \code{P[k,i,j]} is the probability that
#' the system will be in state \code{j} at the \code{k}th time point of
#' \code{t_out} given that it is in state \code{i} at time \code{t_start}.
solve_trans_prob_matrix <- function(system, t_out, log_w0, w = NULL,
                                    log_m = NULL, t_start = 0,
                                    ...) {
  checkmate::assert_class(system, "MultistateSystem")
  checkmate::assert_numeric(t_out)
  K <- length(t_out)
  checkmate::assert_number(t_start, lower = 0)
  H <- system$num_trans()
  W <- system$num_weights()
  S <- system$num_states()
  if (is.null(w)) {
    w <- matrix(0, H, W)
  }
  if (is.null(log_m)) {
    log_m <- rep(0, H)
  }
  kfe <- solve_time_evolution(system, t_out, log_w0, w, log_m, ...)
  P <- array(0, dim = c(K, S, S))
  for (k in seq_len(K)) {
    P[k, , ] <- matrix(kfe[k, 2:ncol(kfe)], S, S)
  }
  P
}

#' Solve transition probability matrices for each subject and draw in
#' 'MultistateModelFit'
#'
#' @param fit A \code{\link{MultistateModelFit}} object
#' @inheritParams solve_trans_prob_matrix
#' @param ... Arguments passed to \code{deSolve::ode()}.
#' @return A list with
#' \itemize{
#'   \item A 4-dimensional array \code{P} where \code{P[n,,,]} is the
#'     \code{P} matrix returned by \code{\link{solve_trans_prob_matrix}} for
#'     subject-draw combination \code{n}
#'    \item Index data frame
#'    \item The numeric vector \code{t_out}
#'  }
solve_trans_prob_matrix_each_subject <- function(fit, t_start = 0, t_out = NULL,
                                                 data = NULL, ...) {
  checkmate::assert_class(fit, "MultistateModelFit")
  checkmate::assert_number(t_start, lower = 0)
  sys <- fit$model$system
  if (is.null(t_out)) {
    t_out <- seq(t_start, sys$get_tmax(), length.out = 30)
  }

  # Get and reshape draws
  S <- fit$num_draws()
  sd <- msmfit_stan_data(fit, data)
  N <- sd$N_sub
  d <- msmfit_inst_hazard_param_draws(fit, data)
  NS <- nrow(d$df)
  pb <- progress::progress_bar$new(total = NS)
  S <- sys$num_states()
  K <- length(t_out)
  A <- array(0, dim = c(NS, K, S, S))
  for (j in seq_len(NS)) {
    pb$tick()
    wj <- d$w[j, , ]
    if (is.null(dim(wj))) {
      wj <- matrix(wj, 1, length(wj))
    }
    A[j, , , ] <- solve_trans_prob_matrix(
      sys, t_out, d$log_w0[j, ], wj, d$log_m[j, ], t_start, ...
    )
  }
  list(P = A, index_df = d$df, t_out = t_out)
}
