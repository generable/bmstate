#' Simulate example data
#'
#' @export
#' @param N number of subjects
#' @param sys_idx index of example system
#' @param h0_base base baseline hazard level
#' @param effect_sizes true effect sizes for the three covariates
#' (\code{age}, \code{sex}, \code{first_dose_amount})
simulate_example_data <- function(N = 10, sys_idx = 1, h0_base = 1e-4,
                                  effect_sizes = c(0.5, 0.5, 0)) {
  possible_covs <- c(
    "sex", "age", "first_dose_amount", "pk_pre_dose", "pk_post_dose",
    "country", "country_num", "dose_adjustment", "dose_arm"
  )
  num_covs <- 3
  NT <- c(4, 9, 25)
  num_trans <- NT[sys_idx]
  num_beta <- num_covs * num_trans
  df_beta_true <- data.frame(
    cov_name = rep(c("age", "sex", "first_dose_amount"), each = num_trans),
    trans_idx = rep(seq_len(num_trans), times = num_covs),
    beta_true = rep(0, num_beta)
  )
  h0_true <- rep(h0_base, num_trans)
  df_beta_true$beta_true[which(df_beta_true$cov_name == "age")] <- effect_sizes[1]
  df_beta_true$beta_true[which(df_beta_true$cov_name == "sex")] <- effect_sizes[2]
  df_beta_true$beta_true[which(df_beta_true$cov_name == "first_dose_amount")] <- effect_sizes[3]
  simulate_multitransition_data(N, possible_covs, h0_true, df_beta_true,
    sys_idx = sys_idx
  )
}


#' Data simulation (multitransition)
#'
#' @param N_subject number of subjects
#' @param covs covariates
#' @param h0_true true baseline hazard value (constant)
#' @param df_beta_true true covariate effects
#' @param sys_idx index of simulation system
#' @param state_names names of the states
simulate_multitransition_data <- function(
    N_subject, covs, h0_true, df_beta_true, sys_idx = 1,
    state_names = NULL) {
  checkmate::assert_integerish(N_subject, lower = 1)
  checkmate::assert_integerish(sys_idx, lower = 1, upper = 3)
  df <- NULL
  pb <- progress::progress_bar$new(total = N_subject)
  n_trans <- length(h0_true)
  m_sub <- matrix(0, N_subject, n_trans)
  for (n in 1:N_subject) {
    pb$tick()
    sex <- n %% 2
    age <- round(30 + 60 * runif(1))
    fda <- c(15, 30, 60)[sample(3, 1)]
    log_hazard_mult <- rep(0, n_trans)
    for (h in 1:n_trans) {
      df_beta_true_h <- df_beta_true |> filter(trans_idx == h)
      log_hazard_mult[h] <- compute_true_hazard_mult(
        sex, age, fda, df_beta_true_h
      )
    }
    m_sub[n, ] <- exp(log_hazard_mult)
    sp <- simulate_path(h0_true, log_hazard_mult, n, sys_idx)
    df_j <- sp$df
    df_j$sex <- sex
    df_j$age <- age
    df_j$country <- 1
    df_j$country_num <- 1
    df_j$first_dose_amount <- fda
    df_j$dose_adjustment <- 0
    df_j$dose_arm <- paste0(df_j$first_dose_amount, ":", df_j$dose_adjustment)
    df_j$pk_pre_dose <- 1
    df_j$pk_post_dose <- 1
    df_j$subject_id <- n
    df <- rbind(df, df_j)
  }

  # Finalize
  term_state <- "Fatality"
  if (is.null(state_names)) {
    if (sys_idx == 1) {
      state_names <- c("Randomization", "Bleed", "Fatality", "Censor")
    } else if (sys_idx == 2) {
      state_names <- c("Randomization", "Bleed", "Stroke", "Fatality", "Censor")
    } else if (sys_idx == 3) {
      state_names <- c(
        "Randomization", "Min Bleed", "Maj Bleed", "Non-Isch Stroke",
        "ICH Bleed", "Fatality", "Censor"
      )
    } else {
      stop("invalid sys_idx")
    }
  }

  list(
    pd = create_sim_pathdata(df, term_state, state_names, covs),
    h0 = sp$conf$h0_all,
    m_sub = m_sub,
    covs = covs,
    df_beta_true = df_beta_true
  )
}

example_3state <- function(t_h0, h0_true, log_hazard_mult) {
  checkmate::assert_numeric(h0_true, len = 4)
  TFI <- matrix(
    c(
      0, 1, 2,
      0, 3, 4,
      0, 0, 0
    ),
    3, 3,
    byrow = TRUE
  )
  h0_ <- rbind(
    t_h0 * 0 + h0_true[1], # rand -> bleed
    t_h0 * 0 + h0_true[2], # bleed -> fatal
    t_h0 * 0 + h0_true[3], # rand -> bleed
    t_h0 * 0 + h0_true[4] # bleed -> fatal
  )
  h0 <- array(0, dim = c(1, nrow(h0_), ncol(h0_)))
  n_trans <- dim(h0)[2]
  h0[1, , ] <- h0_
  m_sub <- array(0, dim = c(1, n_trans, 1))
  m_sub[1, , 1] <- exp(log_hazard_mult)

  # Return
  list(
    TFI = TFI,
    m_sub = m_sub,
    h0 = h0,
    h0_all = h0_true
  )
}

example_4state <- function(t_h0, h0_true, log_hazard_mult) {
  N_trans <- 9
  checkmate::assert_numeric(h0_true, len = N_trans)
  TFI <- matrix(
    c(
      0, 1, 2, 3,
      0, 4, 5, 6,
      0, 7, 8, 9,
      0, 0, 0, 0
    ),
    4, 4,
    byrow = TRUE
  )
  h0_ <- NULL
  for (j in 1:N_trans) {
    h0_ <- rbind(h0_, t_h0 * 0 + h0_true[j])
  }
  h0 <- array(0, dim = c(1, nrow(h0_), ncol(h0_)))
  n_trans <- dim(h0)[2]
  h0[1, , ] <- h0_
  m_sub <- array(0, dim = c(1, n_trans, 1))
  m_sub[1, , 1] <- exp(log_hazard_mult)

  # Return
  list(
    TFI = TFI,
    m_sub = m_sub,
    h0 = h0,
    h0_all = h0_true
  )
}


example_6state <- function(t_h0, h0_true, log_hazard_mult) {
  checkmate::assert_numeric(h0_true, len = 25)
  TFI <- matrix(
    c(
      0, 1, 2, 3, 4, 5,
      0, 6, 7, 8, 9, 10,
      0, 11, 12, 13, 14, 15,
      0, 16, 17, 18, 19, 20,
      0, 21, 22, 23, 24, 25,
      0, 0, 0, 0, 0, 0
    ),
    6, 6,
    byrow = TRUE
  )
  h0_ <- NULL
  for (j in 1:25) {
    h0_ <- rbind(h0_, t_h0 * 0 + h0_true[j])
  }
  h0 <- array(0, dim = c(1, nrow(h0_), ncol(h0_)))
  n_trans <- dim(h0)[2]
  h0[1, , ] <- h0_
  m_sub <- array(0, dim = c(1, n_trans, 1))
  m_sub[1, , 1] <- exp(log_hazard_mult)

  # Return
  list(
    TFI = TFI,
    m_sub = m_sub,
    h0 = h0,
    h0_all = h0_true
  )
}


# Simulate a 3-year path for a single subject
simulate_path <- function(h0_true, log_hazard_mult, subject_idx, sys_idx = 1) {
  init_state <- 1
  t_sub <- 0
  t_start <- 0
  t_max <- 365.25 * 3
  t_h0 <- seq(t_start, t_max, length.out = 300)
  if (sys_idx == 1) {
    conf <- example_3state(t_h0, h0_true, log_hazard_mult)
  } else if (sys_idx == 2) {
    conf <- example_4state(t_h0, h0_true, log_hazard_mult)
  } else if (sys_idx == 3) {
    conf <- example_6state(t_h0, h0_true, log_hazard_mult)
  } else {
    stop("invalid sys_idx")
  }

  h0 <- conf$h0
  TFI <- conf$TFI
  m_sub <- conf$m_sub # hazard multipliers
  n_trans <- dim(h0)[2]
  df <- generate_path(
    init_state = init_state, log_h0 = log(h0),
    t_pred = t_h0, m_sub = log(m_sub), t_sub = t_sub, TFI = TFI,
    t_start = t_start,
    t_max = t_max, n_trans = n_trans,
    discretize = 1
  )
  df <- data.frame(df) |>
    mutate(draw_idx = 1, rep_idx = 1)
  list(
    df = df,
    conf = conf
  )
}

# Hazard multiplier
compute_true_hazard_mult <- function(sex, age, fda, df_beta_true) {
  beta_sex <- df_beta_true |>
    filter(cov_name == "sex") |>
    dplyr::pull(beta_true)
  beta_age <- df_beta_true |>
    filter(cov_name == "age") |>
    dplyr::pull(beta_true)
  beta_dose <- df_beta_true |>
    filter(cov_name == "first_dose_amount") |>
    dplyr::pull(beta_true)
  checkmate::assert_number(beta_sex)
  checkmate::assert_number(beta_age)
  checkmate::assert_number(beta_dose)
  log_hazard_mult <- 0
  if (sex == 1) {
    log_hazard_mult <- log_hazard_mult + beta_sex
  } else {
    log_hazard_mult <- log_hazard_mult - beta_sex
  }
  log_hazard_mult <- log_hazard_mult + beta_age * (age - 60) / 17
  log_hazard_mult <- log_hazard_mult + beta_dose * (fda - 35) / 19
  log_hazard_mult
}


#' Train-test split for a transitions data frame
#'
#' @export
#' @param dt A transitions data frame
#' @return a list with test and train subject indices
do_split <- function(dt) {
  all_sub <- unique(dt$df$subject_id)
  N_sub <- length(all_sub)
  i_test <- sample.int(N_sub, round(N_sub / 4))
  test_sub <- all_sub[i_test]
  train_sub <- setdiff(all_sub, test_sub)
  list(
    test_sub = test_sub,
    train_sub = train_sub
  )
}


# Create final pathdata
create_sim_pathdata <- function(df, term_state, state_names, covs) {
  # Filter censor events after terminal state
  term_state_idx <- which(state_names == term_state)
  df <- df |>
    group_by(subject_id) |>
    mutate(first_termstate_time = ifelse(any(state == term_state_idx),
      min(time[state == term_state_idx]), Inf
    )) |>
    filter(time <= first_termstate_time) |>
    select(-first_termstate_time)

  # Create PathData
  create_pathdata(as_tibble(df), covs,
    state_names = state_names,
    terminal_states = term_state
  )
}
