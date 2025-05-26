#' Solve the Kolmogorov forward equation of a Markovian multistate system
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param t A time grid
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_w0 A vector of length \code{n_trans}
#' @param log_m A vector of length \code{n_trans}
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
#' @param t_init Initial time
#' @param t_end End time
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_m A vector of length \code{n_trans}
#' @return A matrix \code{P} where \code{P[i,j]} is the probability that
#' the system will be in state \code{j} at time \code{t_end} given that it
#' is in state \code{i} at time \code{t_init}
solve_trans_prob_matrix <- function(system, log_w0, w = NULL,
                                    log_m = NULL, t_init = 0, t_end = NULL) {
  checkmate::assert_class(system, "MultistateSystem")
  if (is.null(t_end)) {
    t_end <- system$get_tmax()
  }
  checkmate::assert_number(t_init, lower = 0)
  checkmate::assert_number(t_end, lower = t_init + 1e-9)
  H <- system$num_trans()
  W <- system$num_weights()
  S <- system$num_states()
  if (is.null(w)) {
    w <- matrix(0, H, W)
  }
  if (is.null(log_m)) {
    log_m <- rep(0, H)
  }
  t <- c(t_init, t_end)
  kfe <- solve_time_evolution(system, t, log_w0, w, log_m)
  P <- matrix(kfe[2, 2:ncol(kfe)], S, S)
  cn <- system$tm()$states
  colnames(P) <- cn
  rownames(P) <- cn
  P
}
