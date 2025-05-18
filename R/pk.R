#' Partially steady-state PK model
#'
#' @export
#' @description For each subject
#' @param t vector of output time points
#' @param dose_ss dose amount in SS (for each subject)
#' @param times time points after \code{t_last_ss}
#' @param doses doses taken after \code{t_last_ss}
#' @param theta PK params for each subject
#' @param tau Dosing interval (same for all subjects).
pop_2cpt_partly_ss <- function(t, dose_ss, times, doses, theta, tau) {
  ensure_exposed_stan_functions()
  checkmate::assert_number(tau, lower = 0)
  checkmate::assert_numeric(dose_ss, lower = 0)
  N_sub <- length(dose_ss)
  checkmate::assert_list(times, len = N_sub)
  checkmate::assert_list(doses, len = N_sub)
  checkmate::assert_list(t, len = N_sub)
  checkmate::assert_matrix(theta, nrows = N_sub, ncols = 3)
  N_last <- length(times[[1]])
  amounts <- pop_2cpt_partly_ss_stage1(dose_ss, times, doses, theta, tau)
  pop_2cpt_partly_ss_stage2(t, dose_ss, times, doses, amounts, theta, tau)
}
