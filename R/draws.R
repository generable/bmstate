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
  w <- fit$draws_of("weights") # dim = c(S, H, W)
  log_w0 <- fit$draws_of("log_w0") # dim = c(S, H)
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
