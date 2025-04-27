# Visualize a simple function
plot_lambda <- function(fun, t_min, t_max) {
  lambda_vec <- Vectorize(fun)
  time <- seq(t_min, t_max, length.out = 300)
  rate <- lambda_vec(time)
  df <- data.frame(time, rate)
  ggplot(df, aes(x = time, y = rate)) +
    geom_line() +
    ylim(0, max(rate) * 1.2) +
    xlab("t") +
    ylab("lambda(t)") +
    ggtitle("Rate")
}

# Plot helper
plot_binary_matrix <- function(A, edge = FALSE) {
  df <- melt(A)
  # Rename columns for clarity
  colnames(df) <- c("Row", "Column", "Value")

  # Plot the binary matrix using ggplot2
  plt <- ggplot(df, aes(x = Column, y = Row, fill = factor(Value))) +
    scale_fill_manual(values = c("0" = "gray90", "1" = "gray10")) +
    theme_bw() +
    theme(legend.position = "none")
  if (edge) {
    plt <- plt + geom_tile(color = "white")
  } else {
    plt <- plt + geom_tile() + theme(panel.grid = element_blank())
  }
  plt
}

# Plot single transition
plot_transition <- function(sd, idx) {
  A <- sd$transition[, idx, ]
  plot_binary_matrix(A) + labs(x = "event_idx", y = "subject")
}

# Visualize all transition indices
plot_transitions <- function(sd, names) {
  plots <- list()
  NT <- sd$N_trans
  for (j in seq_len(NT)) {
    plots[[j]] <- plot_transition(sd, j) + ggtitle(names[j])
  }
  ncol <- ceiling(NT / 3)
  ggpubr::ggarrange(plotlist = plots, nrow = 3, ncol = ncol)
}

# Plot used basis functions
plot_bf <- function(sd) {
  t <- rep(sd$t_grid, sd$N_sbf)
  y <- as.vector(sd$SBF_grid)
  idx <- as.factor(rep(1:sd$N_sbf, each = sd$N_grid))
  df <- data.frame(t, y, idx)
  ggplot(df, aes(t, y, color = idx)) +
    geom_line() +
    ggtitle("Basis functions")
}

#' Plot log baseline hazard
#'
#' @inheritParams plot_other_beta
#' @param log_h0_true true baseline hazards
#' @param sd Stan data list
plot_h0 <- function(fit, sd, pd, log_h0_true = NULL) {
  df <- fit |> tidybayes::spread_draws(log_h0[transition, time_idx])
  dt <- pd$as_transitions()
  all_states <- pd$state_names
  if (!is.null(log_h0_true)) {
    df <- df |> left_join(log_h0_true, by = "transition")
  }
  plt <- plot_fun_per_transition(df, sd, legend, "log_h0", "t_pred", all_states) +
    ggtitle("Log baseline hazard") + theme(strip.text = element_text(
      size = 5
    ), legend.position = "none")
  if (!is.null(log_h0_true)) {
    plt <- plt + geom_hline(
      mapping = aes(yintercept = log_h0_true), lty = 3, color = "red", lwd = 1
    )
  }
  plt
}

# Plot age_effect
plot_age_effect <- function(fit, sd, legend, all_states, filter_transitions = NULL,
                            filter_types = NULL) {
  df <- fit |> tidybayes::spread_draws(age_effect[flag, transition, time_idx])
  plot_fun_per_transition(df, sd, legend, "age_effect", "x_age_pred",
    all_states = all_states,
    filter_transitions = filter_transitions, filter_types = filter_types
  ) +
    ggtitle("Age effect")
}

# Plot a function for each transition at t_pred
plot_fun_per_transition <- function(df, sd, legend, name, t_name, all_states,
                                    filter_transitions = NULL, filter_types = NULL) {
  trans_names <- legend$trans_char
  type_names <- .create_types(all_states)
  t_pred <- sd[[t_name]]
  P <- length(t_pred)
  df_t <- data.frame(time_idx = seq_len(P), time = t_pred)
  df <- df |> left_join(df_t, by = "time_idx")
  df$trans_type <- legend$trans_type[df$transition]
  df$Type <- as.factor(type_names[df$trans_type])
  df$Transition <- paste0(df$trans_type, ":", trans_names[df$transition])
  if (!is.null(filter_transitions)) {
    df <- df |>
      filter(Transition %in% filter_transitions)
  }
  if (!is.null(filter_types)) {
    df <- df |>
      filter(Type %in% filter_types)
  }
  df |>
    ggplot(aes(x = time, y = !!sym(name), fill = Type)) +
    stat_lineribbon(alpha = 0.5) +
    facet_wrap(~Transition, ncol = 5) +
    guides(fill = guide_legend(position = "top", nrow = 2))
}

# This can be used after as_transitions
name_legend <- function(legend, state_names) {
  legend$state_char <- state_names[legend$state]
  legend$prev_state_char <- state_names[legend$prev_state]
  s1 <- shorten_name2(legend$prev_state_char)
  s2 <- shorten_name2(legend$state_char)
  legend$trans_char <- paste0(s1, " -> ", s2)
  legend
}

plot_transitions_pd <- function(pd, idx) {
  dat <- pd$df
  all_states <- pd$state_names
  df <- dat |>
    filter(subject_id == unique(dat$subject_id)[idx]) |>
    mutate(
      state_char = all_states[state],
      state = factor(state_char, levels = all_states, ordered = T)
    )
  df$is_event <- as.factor(df$is_event)
  ggplot(df, aes(
    x = time, y = state,
    group = subject_id, pch = is_event
  )) +
    geom_step() +
    geom_point(size = 3, aes(color = is_event)) +
    ggtitle(paste0("Subject ", df$subject_id[1], " (idx = ", idx, ")")) +
    ylab("State") +
    theme_bw()
}

# For debugging Stan data creation for a single subject
debug_standata <- function(pd, dat_trans, sd, idx, trans_names, all_states,
                           id_map_train) {
  checkmate::assert_class(pd, "PathData")
  dat <- pd$as_data_frame()
  sub_id <- subject_idx_to_id(id_map_train, idx)
  df <- dat |>
    filter(subject_id == sub_id) |>
    mutate(
      state_char = all_states[state],
      state = factor(state_char, levels = all_states, ordered = T)
    )
  df$is_event <- as.factor(df$is_event)
  dat_plt <- ggplot(df, aes(
    x = time, y = state,
    group = subject_id, pch = is_event
  )) +
    geom_step() +
    geom_point(size = 3, aes(color = is_event)) +
    ggtitle(paste0("Subject ", df$subject_id[1], " (idx = ", idx, ")")) +
    ylab("State") +
    theme_bw()
  subject_rows <- which(sd$x_sub == idx)
  transition <- sd$transition[, subject_rows, drop = F]
  at_risk <- sd$at_risk[, subject_rows, drop = F]
  rownames(transition) <- trans_names
  rownames(at_risk) <- trans_names
  R <- length(subject_rows)
  t_end <- sd$t_end[subject_rows]
  sd_sub <- list(
    N_sub_int = R,
    transition = transition,
    at_risk = at_risk,
    t_end = t_end # interval end times
  )
  p1 <- plot_binary_matrix(at_risk, edge = TRUE) + ylab("Transition") +
    xlab("Interval") + ggtitle("Risk indicator (at_risk)") +
    scale_x_continuous(breaks = seq_len(R))
  p2 <- plot_binary_matrix(transition, edge = TRUE) + ylab("Transition") +
    xlab("Interval") + ggtitle("Transition indicator (transition)") +
    scale_x_continuous(breaks = seq_len(R))
  mat_plt <- ggarrange(p1, p2, nrow = 1, ncol = 2)
  plt <- ggarrange(dat_plt, mat_plt, nrow = 2, ncol = 1)
  list(
    df = df,
    plt = plt,
    stan_data = sd_sub
  )
}


# Plot KDE of effect multipliers in PK model
plot_effect_beta_pk <- function(a, beta_name, group_by) {
  a[[group_by]] <- as.factor(a[[group_by]])
  ggplot(a, aes(
    x = !!sym(beta_name), group = !!sym(group_by),
    color = !!sym(group_by),
    fill = !!sym(group_by), y = !!sym(group_by)
  )) +
    geom_vline(xintercept = 0, lty = 3) +
    stat_dist_halfeye(alpha = 0.4) +
    theme(legend.position = "none")
}

plot_effects_pk <- function(fit, params) {
  a <- fit |> tidybayes::spread_draws(beta_ka[flag, var])
  a$cov <- params$ka_covariates[a$var]
  p1 <- plot_effect_beta_pk(a, "beta_ka", "cov")
  a <- fit |> tidybayes::spread_draws(beta_CL[flag, var])
  a$cov <- params$CL_covariates[a$var]
  p2 <- plot_effect_beta_pk(a, "beta_CL", "cov")
  a <- fit |> tidybayes::spread_draws(beta_V2[flag, var])
  a$cov <- params$V2_covariates[a$var]
  p3 <- plot_effect_beta_pk(a, "beta_V2", "cov")

  p <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
  annotate_figure(p, top = "Covariate effects in PK model")
}


# Plot KDE of effect multipliers for each transition
# df_true must be a data frame with columns trans_idx and beta_true
plot_effect_beta <- function(a, beta_name, pd, df_true = NULL) {
  dt <- pd$as_transitions()
  legend <- dt$legend
  all_states <- pd$state_names
  type_names <- .create_types(all_states)
  trans_names <- legend$trans_char
  trans_types <- legend$trans_type
  names_long <- paste0(trans_names, " (", 1:length(trans_names), ")")
  if (!is.null(df_true)) {
    a <- a |> left_join(df_true, by = "trans_idx")
  }
  a$Transition <- names_long[a$trans_idx]
  a$trans_type <- trans_types[a$trans_idx]
  a$Type <- as.factor(type_names[a$trans_type])
  a <- a |> group_by(Type)
  plt <- a |>
    ungroup() |>
    ggplot(aes(
      x = !!sym(beta_name), group = trans_idx, color = Type, fill = Type,
      y = Transition
    )) +
    facet_grid(Type ~ ., scale = "free_y", space = "free", switch = "y") +
    ylab("") +
    geom_vline(xintercept = 0, colour = "grey20", linetype = "solid") +
    stat_dist_halfeye(alpha = 0.5) +
    theme(legend.position = "none")
  if (!is.null(df_true)) {
    plt <- plt + geom_vline(
      mapping = aes(xintercept = beta_true), lty = 2,
      color = "steelblue3"
    )
  }

  plt
}

.create_types <- function(states) {
  type_names <- states[2:length(states)]
}

# Plot KDE of effect multipliers for each transition
plot_effect_beta_group <- function(a, var_name, beta_name, legend, all_states) {
  trans_names <- legend$trans_char
  trans_types <- legend$trans_type
  type_names <- .create_types(all_states)
  names_long <- paste0(trans_names, " (", 1:length(trans_names), ")")
  a$Transition <- names_long[a$trans_idx]
  a$trans_type <- trans_types[a$trans_idx]

  a$grp <- paste0(a$trans_idx, a[["var_name"]])
  a$Type <- as.factor(type_names[a$trans_type])
  a |>
    group_by(Type) |>
    mutate(type_median = median(!!sym(beta_name))) |>
    ungroup() |>
    ggplot(aes(
      x = !!sym(beta_name), group = grp, color = Type, fill = Type,
      y = Transition
    )) +
    facet_grid(
      rows = vars(Type), cols = vars(!!sym(var_name)), scale = "free_y",
      space = "free", switch = "y"
    ) +
    ylab("") +
    geom_vline(xintercept = 0, colour = "grey20", linetype = "solid") +
    stat_dist_halfeye(alpha = 0.5) +
    geom_vline(aes(xintercept = type_median, color = Type, group = NULL),
      linetype = "dotted"
    ) +
    theme(legend.position = "none")
}


# Visualize instant hazards for a given hazard multiplier
plot_inst_haz <- function(log_C, h0, sd, trans_names, draw_idx_hl = 1) {
  S <- dim(log_C)[1]
  J <- dim(log_C)[2]
  df <- NULL
  for (j in seq_len(J)) {
    haz <- h0[, j, ] * exp(log_C[, j, 1])
    df_j <- data.frame(
      inst_haz = as.vector(haz),
      time = rep(sd$t_pred, each = S)
    )
    df_j$draw_idx <- rep(1:S, length(sd$t_pred))
    df_j$transition <- j
    df <- rbind(df, df_j)
  }
  df$trans_name <- as.factor(trans_names[df$transition])
  df$hl <- df$draw_idx == draw_idx_hl
  ttl <- paste0("Instant hazards")
  ggplot(df, aes(
    x = time, y = inst_haz, group = draw_idx
  )) +
    geom_line(alpha = 0.25) +
    geom_line(data = df |> filter(hl), color = "red") +
    ylab("Hazard") +
    facet_wrap(. ~ trans_name) +
    scale_y_log10() +
    ggtitle(paste0("highlighting draw ", draw_idx_hl))
}


# PK model fit check
plot_pred_conc <- function(fit, conc, conc_lasttwo, sd, sd_gq, sub_idx,
                           oos = FALSE, cut = TRUE) {
  txt <- "In-sample"
  if (oos) {
    txt <- "Out-of-sample"
  }
  sd_gq <- update_stan_data_pk_pred(sd)
  if (oos) {
    ss_trough <- median(rv(fit, "ss_trough_oos")[sub_idx])
    ss_peak <- median(rv(fit, "ss_peak_oos")[sub_idx])
    ss_auc <- median(rv(fit, "ss_auc_oos")[sub_idx])
    theta_pk <- round(median(rv(fit, "theta_pk_oos"))[1, sub_idx, ], 3)
  } else {
    ss_trough <- median(rv(fit, "ss_trough")[sub_idx])
    ss_peak <- median(rv(fit, "ss_peak")[sub_idx])
    ss_auc <- median(rv(fit, "ss_auc")[sub_idx])
    theta_pk <- round(median(rv(fit, "theta_pk"))[1, sub_idx, ], 3)
  }
  sub <- paste0("median theta = [", paste(theta_pk, collapse = ", "), "]")
  ci <- conc[sub_idx, ]


  if (oos) {
    t <- sd_gq$t_pred_pk_oos[sub_idx, ]
    t_obs <- sd$t_obs_pk_oos[sub_idx, ]
    conc_obs <- sd$conc_pk_oos[sub_idx, ]
    last_t <- sd$last_times_oos[sub_idx, ]
  } else {
    t <- sd_gq$t_pred_pk[sub_idx, ]
    t_obs <- sd$t_obs_pk[sub_idx, ]
    conc_obs <- sd$conc_pk[sub_idx, ]
    last_t <- sd$last_times[sub_idx, ]
    ci_l2 <- conc_lasttwo[sub_idx, ]
    t_ss_end_l2 <- sd$last_two_times[sub_idx, 1]
  }

  df <- data.frame(
    t = t,
    conc = as.vector(median(ci)),
    lower = as.vector(quantile(ci, probs = 0.05)),
    upper = as.vector(quantile(ci, probs = 0.95))
  )
  compare <- !oos
  if (compare) {
    df_l2 <- data.frame(
      t = t,
      conc = as.vector(median(ci_l2)),
      lower = as.vector(quantile(ci_l2, probs = 0.05)),
      upper = as.vector(quantile(ci_l2, probs = 0.95))
    )
    if (cut) {
      df_l2 <- df_l2 |> filter(t >= t_ss_end_l2)
    }
  }
  df_dat <- data.frame(t = t_obs, conc = conc_obs)
  plt <- ggplot(df, aes(x = t, y = conc, ymin = lower, ymax = upper)) +
    geom_vline(
      xintercept = last_t, lty = 2, color = "firebrick"
    ) +
    geom_hline(
      yintercept = c(ss_trough, ss_peak, ss_auc / 24),
      lty = 2, color = "steelblue"
    )
  if (compare) {
    plt <- plt + geom_ribbon(
      data = df_l2,
      alpha = 0.7, fill = "firebrick"
    ) +
      geom_line(data = df_l2, color = "firebrick4", lwd = 1)
  }
  plt <- plt +
    geom_ribbon(alpha = 0.7, fill = "steelblue") +
    geom_line(color = "steelblue4") +
    geom_point(
      data = df_dat, color = "red", mapping = aes(x = t, y = conc),
      inherit.aes = FALSE
    ) +
    ylab("Concentration") +
    xlab("Time (hours)") +
    ggtitle(paste0("Subject ", sub_idx, " (", txt, ")"),
      subtitle = sub
    )
  plt
}

# Annotated plot of predicted concentration for many subjects
plot_pred_conc_many <- function(fit, gq, sd, sd_gq, ids_plot, id_map,
                                oos = FALSE,
                                nrow = 2, ncol = 2, cut = TRUE) {
  plt_conc <- list()
  conc <- rv(gq, "conc_mu_pk_pred")
  conc_l2 <- rv(gq, "conc_mu_pk_pred_lasttwo")
  if (oos) {
    conc <- rv(gq, "conc_mu_pk_pred_oos")
    conc_l2 <- NULL
  }
  j <- 0
  for (is in ids_plot) {
    j <- j + 1
    plt_conc[[j]] <- plot_pred_conc(fit, conc, conc_l2, sd, sd_gq, is, oos, cut) +
      labs(caption = paste0("subject_id = ", subject_idx_to_id(id_map, is)))
  }
  mu_pk <- median(rv(fit, "log_mu_pk"))
  sig_pk <- median(rv(fit, "log_sig_pk"))
  st1 <- paste0("med. log_mu_pk = [", paste(round(mu_pk, 3), collapse = ","), "]")
  st2 <- paste0("med. log_sig_pk = [", paste(round(sig_pk, 3), collapse = ","), "]")
  supertitle <- paste(st1, st2, sep = ", ")
  viz_conc <- ggarrange(plotlist = plt_conc, nrow = nrow, ncol = ncol)
  annotate_figure(viz_conc, top = supertitle)
}

# Correlation plot
plot_cor <- function(sd, x, y, name_x, name_y) {
  df <- data.frame(x, y)
  df$ss_dose <- as.factor(sd$dose_ss)

  # Calculate the correlation coefficient
  correlation <- cor(df$x, df$y)

  # Create the ggplot
  p <- ggplot(df, aes(x = x, y = y, color = ss_dose)) +
    geom_smooth(method = "lm", se = FALSE, color = "gray20", lty = 1) +
    geom_point(alpha = 0.5) +
    annotate("text",
      x = max(df$x), y = max(df$y),
      label = paste("correlation =", round(correlation, 2)),
      hjust = 1.5, vjust = 1, size = 4, color = "gray20"
    ) + # Add text annotation for correlation coefficient
    labs(
      title = paste0(name_x, " vs. ", name_y),
      x = name_x,
      y = name_y
    ) +
    theme_minimal() +
    theme(legend.position = "top")
  p
}

#' Plot coefficients for other covariates
#'
#' @export
#' @param fit fit object
#' @param stan_dat stan data list
#' @param pd A \code{\link{PathData}} object used to create the Stan data.
#' @param df_beta_true must be a data frame with columns \code{trans_idx},
#' \code{cov_name} and \code{beta_true}
plot_other_beta <- function(fit, stan_dat, pd,
                            df_beta_true = NULL) {
  checkmate::assert_class(pd, "PathData")
  names <- stan_dat$x_oth_names
  a <- fit |> tidybayes::spread_draws(beta_oth[cov_idx, trans_idx])
  a$cov_name <- as.factor(names[a$cov_idx])
  plt_oth <- NULL
  for (j in seq_len(stan_dat$stan_data$N_oth)) {
    a_nam <- a |> filter(cov_name == names[j])
    if (!is.null(df_beta_true)) {
      df_true <- df_beta_true |> filter(cov_name == names[j])
    } else {
      df_true <- NULL
    }
    plt_oth[[j]] <- plot_effect_beta(
      a_nam, "beta_oth", pd, df_true
    ) + ggtitle(names[j]) + xlab(expression(beta))
  }
  plt_oth
}

#' Plot event rates over all subjects
#'
#' @export
#' @param pd A \code{\link{PathData}} object of observed data
#' @param paths_gen A \code{\link{PathData}} object of generated paths
#' (in-sample)
#' @param paths_gen_oos A \code{\link{PathData}} object of generated paths
#' (out-of-sample)
#' @param names names of the three former objects
plot_p_event <- function(pd, paths_gen, paths_gen_oos, names) {
  df1 <- p_event(pd$path_df, pd$state_names) # Observed
  df2 <- p_event(paths_gen$path_df, paths_gen$state_names) # Pred.
  df3 <- p_event(paths_gen_oos$path_df, paths_gen_oos$state_names) # Pred.(oos)

  df1$name <- names[1]
  df2$name <- names[2]
  df3$name <- names[3]
  df <- rbind(df1, df2, df3)
  df$name <- as.factor(df$name)
  ggplot(
    df,
    aes(y = p_paths, x = state_char, group = name, fill = name)
  ) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab("Event") +
    ylab("Proportion of paths") +
    ggtitle("Proportion of paths with given event") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top", legend.title = element_blank()
    ) +
    scale_fill_brewer(palette = 3, type = "qual")
}

#' Plot event rates over all subjects by dose
#'
#' @export
#' @inheritParams plot_p_event
#' @param t_max max time used in path generation
plot_p_event_by_dose <- function(pd, paths_gen, paths_gen_oos, names,
                                 t_max = 365.25 * 3) {
  pe_data <- summarize_ppsurv(pd,
    by = c("first_dose_amount"),
    target_times = t_max
  ) |>
    mutate(
      p_paths = 1 - surv,
      state_char = pd$state_names[state]
    )
  pe_pred <- summarize_ppsurv(paths_3yr,
    by = c("first_dose_amount"),
    target_times = t_max
  ) |>
    mutate(
      p_paths = 1 - surv,
      state_char = pd$state_names[state]
    )
  pe_oos <- summarize_ppsurv(paths_3yr_oos,
    by = c("first_dose_amount"),
    target_times = t_max
  ) |>
    mutate(
      p_paths = 1 - surv,
      state_char = pd$state_names[state]
    )
  pe_data$name <- names[1]
  pe_pred$name <- names[2]
  pe_oos$name <- names[3]

  res <- rbind(pe_data, pe_pred, pe_oos)
  res$state_char <- as.factor(res$state_char)
  ggplot(res, aes(x = first_dose_amount, y = p_paths, color = state_char)) +
    geom_line() +
    geom_point() +
    ylab("Event probability") +
    xlab("Assigned dose (mg)") +
    theme(legend.title = element_blank()) +
    facet_wrap(. ~ name) +
    theme(legend.position = "top")
}

# Visualize PK model residuals
create_pk_residual_df <- function(fit, sd, oos = FALSE, pre = TRUE) {
  if (pre) {
    idx <- 1
    type <- "Pre-dose"
  } else {
    idx <- 2
    type <- "Post-dose"
  }
  if (!oos) {
    ts <- sd$t_since_last_pk[, idx]
    nam <- "conc_mu_pk"
    ynam <- "conc_pk"
    mode <- "In-sample"
  } else {
    ts <- sd$t_since_last_pk_oos[, idx]
    nam <- "conc_mu_pk_oos"
    ynam <- "conc_pk_oos"

    mode <- "Out-of-sample"
  }
  mu <- mean(rv(fit, nam))[1, , idx]
  y <- sd[[ynam]][, idx]
  residual <- y - mu
  df <- data.frame(t_since_dose = ts, residual = residual)
  df$mode <- mode
  df$type <- type
  df
}

# Visualize PK model residuals
plot_pk_residual <- function(fit, sd) {
  df <- rbind(
    create_pk_residual_df(fit, sd, oos = FALSE, pre = TRUE),
    create_pk_residual_df(fit, sd, oos = TRUE, pre = TRUE),
    create_pk_residual_df(fit, sd, oos = FALSE, pre = FALSE),
    create_pk_residual_df(fit, sd, oos = TRUE, pre = FALSE)
  )
  ggplot(as_tibble(df), aes(x = t_since_dose, y = residual, color = mode)) +
    geom_point(alpha = 0.75, pch = 4) +
    facet_wrap(. ~ type) +
    theme(legend.position = "top")
}

# Plot instant hazards for one subject
plot_inst_haz_subject <- function(fit, sd, sub_idx, draw_idx, oos, tr_names) {
  checkmate::assert_number(draw_idx)
  checkmate::assert_number(sub_idx)
  h0 <- exp(get_and_format_log_h0_draws(fit))
  i1_idx <- sd$sub_start_idx[sub_idx]
  log_C <- get_and_format_log_C_haz_draws(fit, oos)[, , i1_idx, drop = F]
  plot_inst_haz(log_C, h0, sd, tr_names, draw_idx) +
    ggtitle(paste0("subject = ", sub_idx, ", draw = ", draw_idx))
}

# Debug
plot_inst_haz_longest_path <- function(fit, sd, paths, id_map, oos, tr_names) {
  lp <- paths$longest_path()
  df <- lp$as_data_frame()
  len <- nrow(df)
  message("Longest path has ", len, " transitions")
  sid <- unique(df$subject_id)[1] # should have length 1
  draw_idx <- as.numeric(unique(df$draw_idx))[1] # should have length 1
  sub_idx <- which(id_map$subject_id == sid)
  plot_inst_haz_subject(fit, sd, sub_idx, draw_idx, oos, tr_names)
}

#' Plot Brier scores
#'
#' @export
#' @param ppsurv_subj \code{ppsurv} object
#' @param pd \code{\link{PathData}} object
#' @param sub_ids_char_training training subject ids (character)
create_brier_score_plot <- function(ppsurv_subj, pd, sub_ids_char_training) {
  checkmate::assert_class(pd, "PathData")
  sub_ids_char_test <- unique(ppsurv_subj$subject_id)
  scores <- compute_scores(
    pd,
    ppsurv = ppsurv_subj,
    which = "brier_score",
    keep_gt_info = FALSE,
    sub_ids_char = sub_ids_char_test, # set by default from ppsurv_subj
    sub_ids_char_training = sub_ids_char_training
  )

  # no need to recompute score for predicted paths
  km_scores_vs_dose <- compute_scores(
    pd,
    target_times = sort(unique(ppsurv_subj$time)),
    km_formula = ~dose_arm,
    which = "brier_score",
    keep_gt_info = FALSE,
    sub_ids_char = sub_ids_char_test, # needs explicit test set
    sub_ids_char_training = sub_ids_char_training
  )

  # plot
  scores |>
    mutate(method = if_else(method == "km", "km overall", method)) |>
    bind_rows(km_scores_vs_dose |> mutate(method = "km per dose")) |>
    mutate(
      Event = factor(str_wrap(pd$state_names[state], width = 10, whitespace_only = T)),
      Event = forcats::fct_reorder(Event, state)
    ) |>
    ggplot(aes(
      x = time / 365.25, y = brier_score,
      group = str_c(Event, method),
      colour = method, fill = method
    )) +
    geom_line() +
    facet_wrap(~Event, scale = "free_y") +
    scale_y_continuous("Brier Score") +
    scale_x_continuous("Study Time (years)")
}

#' Plot concordance indices
#'
#' @export
#' @param ppsurv_subj \code{ppsurv} object
#' @param ppsurv_subj_oos \code{ppsurv} object (out-of-sample)
#' @param sub_ids_train training subject ids
#' @param sub_ids_test test subject ids
#' @inheritParams plot_p_event
create_cindex_plot <- function(pd, paths_gen, paths_gen_oos,
                               sub_ids_train, sub_ids_test,
                               ppsurv_subj, ppsurv_subj_oos) {
  all_states <- pd$state_names

  # Implementation 1
  ci_is <- list()
  ci_oos <- list()
  ci_is_table <- NULL
  ci_oos_table <- NULL
  km_fits <- list()
  t_eval <- 3 * 365.25
  for (j in seq_len(length(all_states) - 2)) {
    state_idx <- j + 1
    km_fits[[j]] <- km_fit(pd, state_idx, sub_ids_train)
    ci_is[[j]] <- c_index(
      pd, paths_3yr, state_idx, km_fits[[j]], sub_ids_train, t_eval
    )
    ci_oos[[j]] <- c_index(
      pd, paths_3yr_oos, state_idx, km_fits[[j]], sub_ids_test, t_eval
    )
    ci_is_table <- rbind(ci_is_table, ci_is[[j]]$table)
    ci_oos_table <- rbind(ci_oos_table, ci_oos[[j]]$table)
  }
  names(ci_is) <- sapply(ci_is, function(x) x$state_name)
  names(ci_oos) <- sapply(ci_oos, function(x) x$state_name)

  # Implementation 2
  ci_is2 <- compute_scores(
    pd,
    ppsurv = ppsurv_subj |> filter(time == 365.25 * 3),
    km_formula = ~dose_arm,
    which = "cindex",
    keep_gt_info = FALSE
  )

  ci_oos2 <- compute_scores(
    pd,
    ppsurv = ppsurv_subj_oos |> filter(time == 365.25 * 3),
    sub_ids_char_training = sub_ids_train,
    km_formula = ~dose_arm,
    which = "cindex",
    keep_gt_info = FALSE
  )

  # Comparison table
  ci_is_table <- ci_is_table |>
    left_join(
      ci_is2 |> select(method, Event, cindex),
      by = join_by(method, Event),
      suffix = c(".v1", ".v2")
    ) |>
    select(method, Event, cindex.v1, cindex.v2)
  ci_oos_table <- ci_oos_table |>
    left_join(
      ci_oos2 |> select(method, Event, cindex),
      by = join_by(method, Event),
      suffix = c(".v1", ".v2")
    ) |>
    select(method, Event, cindex.v1, cindex.v2)

  # Plot
  plt_ci_is <- annotate_figure(plot_ci_multi(ci_is), top = "In-sample")
  plt_ci_oos <- annotate_figure(plot_ci_multi(ci_oos), top = "Out-of-sample")
  dplyr::lst(ci_is_table, ci_oos_table, plt_ci_is, plt_ci_oos)
}
