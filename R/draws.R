# Utility
rv <- function(fit, name) {
  posterior::as_draws_rvars(fit$draws(name))[[name]]
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


# Get log hazard multiplier of each subject
msmsf_log_m_per_subject <- function(fit) {
  log_m <- fit$draws_of("log_C_haz")
  idx_sub <- fit$stan_data$idx_sub
  N_sub <- fit$stan_data$N_sub
  first_indices <- sapply(seq_len(N_sub), function(x) which(idx_sub == x)[1])
  log_m[, first_indices, , drop = FALSE]
}

#' Extract and reshape draws of instant hazard related parameters
#'
#' @export
#' @inheritParams generate_paths
#' @return a list with elements \code{log_m}, \code{log_w0}, \code{w}, each
#' of which is an array where the first dimension is number of subjects times
#' number of draws
get_inst_hazard_param_draws <- function(fit) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  S <- fit$num_draws()
  N <- fit$stan_data$N_sub
  w <- fit$draws_of("weights")
  log_w0 <- fit$draws_of("log_w0")
  w_rep <- abind::abind(replicate(N, w, simplify = FALSE), along = 1)
  log_w0_rep <- abind::abind(replicate(N, log_w0, simplify = FALSE), along = 1)
  log_m <- msmsf_log_m_per_subject(fit)
  H <- dim(log_m)[3]
  log_m_reshaped <- matrix(aperm(log_m, c(1, 2, 3)), nrow = N * S, ncol = H)
  list(
    log_m = log_m_reshaped,
    log_w0 = log_w0_rep,
    w = w_rep,
    subject_index = rep(seq_len(N), times = S),
    draw_index = rep(seq_len(S), each = N)
  )
}
