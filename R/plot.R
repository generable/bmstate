# Plot KDE of effect multipliers in PK model
plot_effect_beta_pk <- function(a, beta_name, group_by) {
  a[[group_by]] <- as.factor(a[[group_by]])
  ggplot(a, aes(
    x = !!sym(beta_name), group = !!sym(group_by),
    color = !!sym(group_by),
    fill = !!sym(group_by), y = !!sym(group_by)
  )) +
    geom_vline(xintercept = 0, lty = 3) +
    stat_halfeye(alpha = 0.4) +
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

  p <- ggpubr::ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
  ggpubr::annotate_figure(p, top = "Covariate effects in PK model")
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
  a <- a |> dplyr::group_by(Type)
  plt <- a |>
    dplyr::ungroup() |>
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

# PK model fit check
plot_pred_conc <- function(fit, conc, conc_lasttwo, sd, sd_gq, sub_idx,
                           oos = FALSE, cut = TRUE) {
  txt <- "In-sample"
  if (oos) {
    txt <- "Out-of-sample"
  }
  sd_gq <- update_stan_data_pk_pred(sd)
  if (oos) {
    ss_trough <- stats::median(rv(fit, "ss_trough_oos")[sub_idx])
    ss_peak <- stats::median(rv(fit, "ss_peak_oos")[sub_idx])
    ss_auc <- stats::median(rv(fit, "ss_auc_oos")[sub_idx])
    theta_pk <- round(stats::median(rv(fit, "theta_pk_oos"))[1, sub_idx, ], 3)
  } else {
    ss_trough <- stats::median(rv(fit, "ss_trough")[sub_idx])
    ss_peak <- stats::median(rv(fit, "ss_peak")[sub_idx])
    ss_auc <- stats::median(rv(fit, "ss_auc")[sub_idx])
    theta_pk <- round(stats::median(rv(fit, "theta_pk"))[1, sub_idx, ], 3)
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
    conc = as.vector(stats::median(ci)),
    lower = as.vector(stats::quantile(ci, probs = 0.05)),
    upper = as.vector(stats::quantile(ci, probs = 0.95))
  )
  compare <- !oos
  if (compare) {
    df_l2 <- data.frame(
      t = t,
      conc = as.vector(stats::median(ci_l2)),
      lower = as.vector(stats::quantile(ci_l2, probs = 0.05)),
      upper = as.vector(stats::quantile(ci_l2, probs = 0.95))
    )
    if (cut) {
      df_l2 <- df_l2 |> dplyr::filter(t >= t_ss_end_l2)
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
  mu_pk <- stats::median(rv(fit, "log_mu_pk"))
  sig_pk <- stats::median(rv(fit, "log_sig_pk"))
  st1 <- paste0("med. log_mu_pk = [", paste(round(mu_pk, 3), collapse = ","), "]")
  st2 <- paste0("med. log_sig_pk = [", paste(round(sig_pk, 3), collapse = ","), "]")
  supertitle <- paste(st1, st2, sep = ", ")
  viz_conc <- ggpubr::ggarrange(plotlist = plt_conc, nrow = nrow, ncol = ncol)
  ggpubr::annotate_figure(viz_conc, top = supertitle)
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
    a_nam <- a |> dplyr::filter(cov_name == names[j])
    if (!is.null(df_beta_true)) {
      df_true <- df_beta_true |> dplyr::filter(cov_name == names[j])
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
  pe_pred <- summarize_ppsurv(paths_gen,
    by = c("first_dose_amount"),
    target_times = t_max
  ) |>
    mutate(
      p_paths = 1 - surv,
      state_char = pd$state_names[state]
    )
  pe_oos <- summarize_ppsurv(paths_gen_oos,
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


# Visualize a graph
transition_matrix_plot <- function(f, terminal_states, null_state,
                                   edge_labs, col_term = "firebrick",
                                   col_null = "steelblue2", ...) {
  cn <- colnames(f)
  idx_term <- which(colnames(f) %in% terminal_states)
  idx_null <- which(colnames(f) %in% null_state)
  N <- length(cn)
  color <- rep("black", N)
  color[idx_term] <- col_term
  color[idx_null] <- col_null
  acol <- matrix("black", nrow(f), ncol(f))
  lcol <- acol
  lcol[, idx_term] <- col_term
  lcol[, idx_null] <- col_null

  # Create plot
  qgraph::qgraph(f,
    edge.labels = edge_labs, label.color = color,
    edge.color = acol,
    fade = FALSE,
    layout = "circle",
    ...
  )
}
