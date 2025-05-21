# Utility
rv <- function(fit, name) {
  posterior::as_draws_rvars(fit$draws(name))[[name]]
}

#' Extract log hazard multiplier draws
#'
#' @export
#' @param fit fit object
#' @param oos for out-of-sample subjects?
#' @return an array with \code{dim = c(n_draws, n_hazards, n_obs)}
get_and_format_log_C_haz_draws <- function(fit, oos) {
  vn <- "log_C_haz"
  if (oos) {
    vn <- paste0(vn, "_oos")
  } else {
    vn <- paste0(vn, "_is")
  }
  log_C <- rv(fit, vn)
  log_C <- posterior::merge_chains(log_C) # shape c(N_obs, N_pred)
  N_draws <- posterior::ndraws(log_C)
  N_obs <- dim(log_C)[1]
  H <- dim(log_C)[2]
  out <- array(0, dim = c(N_draws, H, N_obs))
  for (j in seq_len(H)) {
    log_C_j <- posterior::as_draws_array(as.vector(log_C[, j]))
    for (d in 1:N_draws) {
      out[d, j, ] <- as.vector(log_C_j[d, 1, ])
    }
  }
  out
}
