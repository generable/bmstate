#' Visualize parts of the 'Stan' data
#'
#' @export
#' @param model A \code{\link{MultistateModel}}.
#' @param stan_data A list of data for Stan. Created using
#' \code{\link{create_stan_data}}.
#' @param name Name of Stan data field to plot. Currently must be one of
#' \code{"transition", "at_risk"}.
#' @param subject_idx Index of subject whose data rows to plot.
plot_stan_data <- function(model, stan_data, name, subject_idx) {
  checkmate::assert_integerish(subject_idx, len = 1)
  checkmate::assert_class(model, "MultistateModel")
  opts <- c("transition", "at_risk")
  checkmate::assert_character(name, len = 1, min.chars = 2)
  checkmate::assert_choice(name, opts)
  plot_intervals_mat(model, stan_data, name, subject_idx)
}

plot_intervals_mat <- function(model, sd, name, subject_idx) {
  rows <- which(sd$idx_sub == subject_idx)
  if (length(rows) == 0) {
    stop("subject not found")
  }
  names <- model$system$tm()$trans_df()$trans_char
  mat <- sd[[name]]
  mat <- mat[rows, , drop = FALSE]
  K <- ncol(mat)
  mat <- mat[, 1:(K - 1)]
  rownames(mat) <- names
  colnames(mat) <- paste(": ", 1:(K - 1))
  mat <- reshape2::melt(mat)
  colnames(mat) <- c("Transition", "Interval", "Value")
  plot_mat(mat) + ggtitle(name, subtitle = "black = 1, white = 0")
}

plot_mat <- function(x) {
  ggplot(x, aes(.data$Interval, .data$Transition, alpha = .data$Value)) +
    geom_tile(colour = "gray50") +
    scale_alpha_identity(guide = "none") +
    coord_equal(expand = 0) +
    theme_bw() +
    theme(panel.grid.major = element_blank())
}

plot_integral <- function(mod, sd, iidx, h0, w) {
  checkmate::assert_number(h0, lower = 0)
  if (is.null(w)) {
    w <- stats::rnorm(mod$system$num_weights(), sd = 0.3)
  }
  checkmate::assert_vector(w)
  hfun <- function(t) {
    log_h <- mod$system$log_baseline_hazard(t, log(h0), w = w)
    exp(log_h)
  }
  tt <- seq(0, mod$get_tmax())
  ll <- hfun(tt)
  t_eval <- sd$t_grid
  l_eval <- hfun(t_eval)
  i1 <- sd$t_start_idx_m1[iidx]
  i2 <- sd$t_end_idx[iidx] - 1
  xint <- c(sd$t_int_start[iidx], sd$t_int_end[iidx])
  ts <- paste0("interval_idx = ", iidx, ", correction_multiplier = ", round(sd$correction_multiplier[iidx], 5))
  ggplot(data.frame(time = tt, hazard = ll), aes(x = time, y = hazard)) +
    stat_function(
      fun = hfun,
      xlim = xint,
      geom = "area",
      fill = "steelblue"
    ) +
    geom_col(
      data = data.frame(time = t_eval, hazard = l_eval),
      fill = NA, col = "gray40", width = delta_grid
    ) +
    geom_col(
      data = data.frame(time = t_eval[i1:i2], hazard = l_eval[i1:i2]),
      fill = "firebrick3", col = "firebrick3", width = delta_grid, alpha = 0.3
    ) +
    geom_line(lwd = 1) +
    ggtitle(ts, subtitle = "Blue area ~ Area of red bars x Correction multiplier")
}
