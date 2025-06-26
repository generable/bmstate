# Options
params <- list(
  N_subject = 20,
  iter_warmup = 400,
  iter_sampling = 200,
  chains = 1,
  n_rep = 100
)

test_that("entire workflow works with more complex model", {
  # Create transition matrix
  mat <- matrix(0, 4, 4)
  mat[1, 2] <- 1
  mat[2, 3] <- 1
  mat[3, 2] <- 1
  mat[3, 4] <- 1
  mat[2, 4] <- 1
  mat[1, 4] <- 1
  mat[1, 3] <- 1
  mat[2, 2] <- 0
  mat[3, 3] <- 0
  sn <- c("Randomization", "Bleed", "Stroke", "Death")
  tm <- TransitionMatrix$new(mat, sn)
  t3yr <- 3 * 365.25
  mod <- create_msm(tm,
    hazard_covs = c("age", "dose_amt"), num_knots = 4,
    tmax = t3yr
  )
  tm0 <- transmat_survival(sn[c(1, 4)])
  mod_tte <- create_msm(tm0,
    hazard_covs = c("age", "dose_amt"), num_knots = 4,
    tmax = t3yr
  )

  NTT <- mod$system$num_states() - 1
  bh_true <- matrix(0, NTT, 2)
  bh_true[1, 2] <- -1.5
  bh_true[2, 2] <- 0.5
  bh_true[3, 1] <- 0.3
  sn <- tm$states_df() |>
    dplyr::filter(!source) |>
    dplyr::pull(state)
  rownames(bh_true) <- paste0("Effect on ", sn)
  colnames(bh_true) <- mod$covs()

  h0_true <- 0.5 * 1e-3
  h0_true_vec <- rep(h0_true, 7)
  h0_true_vec[3] <- 0.1 * h0_true
  h0_true_vec[5] <- 5 * h0_true
  h0_true_vec[7] <- 20 * h0_true

  simdat <- mod$simulate_data(params$N_subject,
    beta_haz = bh_true, w0 = h0_true_vec
  )
  df_h0t <- tm$trans_df()
  df_h0t$h0_true <- h0_true_vec
  simdat <- split_data(simdat, p_test = 0.5)
  train_data <- simdat$train
  test_data <- simdat$test
  pth1 <- as_single_event(train_data$paths, "Death")
  pth2 <- as_single_event(test_data$paths, "Death")
  train_data_tte <- JointData$new(pth1, NULL)
  test_data_tte <- JointData$new(pth2, NULL)

  # Fit the model
  fit <- fit_stan(mod, train_data,
    iter_warmup = params$iter_warmup,
    iter_sampling = params$iter_sampling,
    chains = params$chains,
    adapt_delta = 0.95,
    init = 0.1
  )
  fit_tte <- fit_stan(mod_tte, train_data_tte,
    iter_warmup = params$iter_warmup,
    iter_sampling = params$iter_sampling,
    chains = params$chains,
    adapt_delta = 0.95,
    init = 0.1
  )
  fit_pm <- fit$mean_fit()
  fit_tte_pm <- fit_tte$mean_fit()

  plt1 <- fit$plot_h0() + geom_hline(data = df_h0t, aes(yintercept = h0_true), lty = 2)
  plt2 <- fit_pm$plot_h0() + geom_hline(data = df_h0t, aes(yintercept = h0_true), lty = 2)

  df_beta_true <- data.frame(bh_true) |>
    rownames_to_column("event") |>
    pivot_longer(cols = -event, names_to = "covariate", values_to = "beta")

  df_beta <- t(fit$get_draws("beta_oth"))
  rownames(df_beta) <- paste0("Effect on ", sn)
  colnames(df_beta) <- mod$covs()
  df_beta_long <- data.frame(df_beta) |>
    rownames_to_column("event") |>
    pivot_longer(cols = -event, names_to = "covariate", values_to = "beta")
  plt_beta <- ggplot(df_beta_long, aes(y = "", dist = beta)) +
    stat_dist_halfeye(color = "gray10", fill = "firebrick", alpha = 0.8) +
    labs(x = "Coefficient", y = "Posterior density") +
    facet_wrap(. ~ covariate + event) +
    geom_vline(data = df_beta_true, mapping = aes(xintercept = beta), lty = 2)
  plt_beta

  a <- fit_tte$get_draws("beta_oth")
  names(a) <- mod$covs()
  a <- data.frame(t(a))
  rownames(a) <- paste0("Effect on ", "Death")
  a

  predict_death_risk <- function(fit, dat, n_rep = 100) {
    paths <- generate_paths(fit, data = dat, n_rep = n_rep)
    ev <- "Death"
    trans_prob <- solve_trans_prob_fit(fit,
      data = dat
    )[, c("subject_index", ev)]

    er <- event_risk(paths, ev)
    a <- as_survival(test_data$paths, ev)
    a <- a |> dplyr::left_join(er, by = "subject_id")
    a <- a |> dplyr::select(-"path_id")
    a$dose_amt <- as.factor(a$dose_amt)
    a$risk_analytic <- trans_prob[, ev]
    list(
      df = a,
      p = paths,
      tp = trans_prob
    )
  }


  r1 <- predict_death_risk(fit_pm, test_data, params$n_rep)
  r2 <- predict_death_risk(fit_tte_pm, test_data_tte, params$n_rep)

  df1 <- r1$df
  df1$risk <- df1$risk_analytic
  df2 <- r2$df
  df2$risk <- df2$risk_analytic
  a1 <- c_index(df1)
  a2 <- c_index(df2)
  r <- data.frame(
    model = c("Multistate", "Time-to-event"),
    concordance = c(a1$concordance, a2$concordance)
  )
  r

  risk_plots <- function(a1, a2) {
    a <- rbind(a1, a2)
    m1 <- "Multistate model"
    m2 <- "Time-to-event model"
    dlab <- "Predicted death risk"
    a$Model <- c(rep(m1, nrow(a1)), rep(m2, nrow(a2)))
    a$Censored <- as.logical(a$is_censor)
    p1 <- ggplot(a, aes(x = age, y = risk_analytic, color = dose_amt, group = Model)) +
      geom_point(alpha = 0.7) +
      xlab("Age (years)") +
      theme(legend.position = "top") +
      facet_wrap(. ~ Model) +
      ylab(dlab)
    p2 <- ggplot(a, aes(
      x = time, y = risk_analytic, color = dose_amt, group = Model,
      pch = Censored
    )) +
      geom_point(alpha = 0.7) +
      xlab("Age (years)") +
      theme(legend.position = "top") +
      facet_wrap(. ~ Model) +
      xlab("Survival time") +
      ylab(dlab)
    p3 <- ggplot(a, aes(
      x = risk_analytic, y = risk, color = dose_amt, group = Model,
      pch = Censored
    )) +
      geom_point(alpha = 0.7) +
      xlab(dlab) +
      theme(legend.position = "top") +
      facet_wrap(. ~ Model) +
      ylab(paste(dlab, " (paths)")) +
      xlab(dlab) +
      geom_line(
        data = data.frame(x = c(0, 1), y = c(0, 1)), aes(x, y), lty = 2,
        inherit.aes = FALSE
      )
    list(
      a = ggpubr::ggarrange(p1, p2, nrow = 2, labels = "auto"),
      b = p3
    )
  }

  pr <- risk_plots(r1$df, r2$df)
  pr$a

  pr$b

  pf1 <- p_event(r1$p, by = "dose_amt") |>
    filter(state == "Death") |>
    select(dose_amt, prob)
  pf2 <- p_event(r2$p, by = "dose_amt") |> select(dose_amt, prob)
  of <- p_event(test_data$paths, by = "dose_amt") |>
    filter(state == "Death") |>
    select(dose_amt, prob)
  pf <- pf1 |> left_join(pf2, by = c("dose_amt"))
  pf <- pf |> left_join(of, by = c("dose_amt"))
  colnames(pf) <- c("Dose (mg)", "Multistate model", "Time-to-event model", "Observed")
  round(pf, 3)
})
