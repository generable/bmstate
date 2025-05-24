#' Solve the Kolmogorov forward equation of a Markovian multistate system
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param t A time grid
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_w0 A vector of length \code{n_trans}
#' @param log_m A vector of length \code{n_trans}
solve_time_evolution <- function(system, t, w, log_w0, log_m) {
  checkmate::assert_class(system, "MultistateSystem")
  S <- sys$num_trans()
  W <- sys$num_weights()
  checkmate::assert_numeric(t, min.len = 1)
  checkmate::assert_matrix(w, ncols = S, nrows = W)
  checkmate::assert_numeric(log_w0, len = S)
  checkmate::assert_numeric(log_m, len = S)
  H <- system$log_inst_hazard(t, w, log_w0, log_m)
  H
}
