.form_rhs_syms <- function(formula) {
  # confirm that formula contains no response variables
  if (attr(terms(formula), "response") == 1) {
    stop("Provided formula contains a response but should not. Try providing in the format: `~ 1` or `~ variable`")
  }
  # extract RHS variables as symbols
  term_var_syms <- rlang::syms(attr(terms(formula), "variables", exact = TRUE))
  # remove first element, which is `symbol: list`
  return(term_var_syms[-1])
}

.prepare_observed_data <- function(pd, target_times,
                                   sub_ids_char = NULL,
                                   sub_ids_char_test = NULL,
                                   sub_ids_char_training = NULL,
                                   formula = ~1) {
  # parse formulas early to catch input errors
  surv_lhs_formula <- formula(survival::Surv(time, is_event) ~ 1)
  surv_formula <- update(surv_lhs_formula, formula)
  surv_vars <- .form_rhs_syms(formula) # symbols
  surv_covs <- purrr::map_chr(surv_vars, rlang::as_name) # characters


  # base observed events
  observed_events <- pd$as_data_frame(cov = c("subject_id", surv_covs)) |>
    as_time_to_first_event(
      by = c("subject_id", surv_covs),
      states = pd$get_event_states()
    ) |>
    dplyr::select(subject_id, state, time, is_event, !!!surv_vars) |>
    mutate(one = 1)
  if (!is.null(sub_ids_char) && is.null(sub_ids_char_test) &&
    is.null(sub_ids_char_training)) {
    # in this case, test & training sets are the same
    sub_ids_char_test <- sub_ids_char_training <- sub_ids_char
  } else if (
    is.null(sub_ids_char) && !is.null(sub_ids_char_training) &&
      !is.null(sub_ids_char_test)
  ) {
    # in this case, test & training sets are different. the following line is
    # not required
    # including this to simplify the following error message
    sub_ids_char <- sub_ids_char_test
  } else if (!is.null(sub_ids_char)) {
    sub_ids_char_test <- sub_ids_char
  } else {
    stop("You must provide either sub_ids_char or sub_ids_char_test & training, but neither were provided. Pls review your arguments.")
  }
  test_observed_events <- observed_events |>
    dplyr::filter(subject_id %in% sub_ids_char_test)
  training_observed_events <- observed_events |>
    dplyr::filter(subject_id %in% sub_ids_char_training)

  # we always compute true events & censoring distribution using test dataset
  prepared_observed_events <- test_observed_events |>
    dplyr::select(state, time, is_event, subject_id) |>
    dplyr::group_by(state) |>
    tidyr::nest(event_data = c(time, is_event, subject_id)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      event_data = purrr::map(event_data, ~ dplyr::arrange(.x, subject_id))
    ) |>
    dplyr::mutate(
      st_surv = purrr::map(event_data, ~ survival::Surv(.x$time, .x$is_event)),
      gt_times = purrr::map2(st_surv, event_data, ~ .Gt(.x, .y$time))
    )

  # we also prepare similar Pr using KM estimator as a benchmark
  # this should use the training data projected onto the test population,
  # for each state
  .get_surv_strata <- function(km_fit, target_times) {
    ss <- summary(km_fit, target_times)
    if (length(surv_covs) > 0) {
      ss_tbl <- dplyr::as_tibble(ss[c("time", "surv", "strata")]) |>
        rgeco:::.tidy_km_strata()
    } else {
      ss_tbl <- dplyr::as_tibble(ss[c("time", "surv")]) |>
        dplyr::mutate(VAR = 1)
    }
    ss_tbl
  }
  # here we estimate event rate using the training data
  if (length(surv_covs) == 0) {
    surv_vars <- syms("one")
  }
  p_event_groups <- training_observed_events |>
    dplyr::mutate(VAR = dplyr::dense_rank(!!!surv_vars)) |>
    distinct(VAR, !!!surv_vars)
  p_event <- training_observed_events |>
    mutate(one = 1) |>
    left_join(p_event_groups, by = dplyr::join_by(!!!surv_vars)) |>
    dplyr::select(state, subject_id, time, is_event, VAR) |>
    dplyr::group_by(state) |>
    tidyr::nest(event_data = c(subject_id, time, is_event, VAR)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      # Surv estimator per state
      event_surv = furrr::future_map(
        event_data,
        ~ survival::survfit(update(surv_lhs_formula, ~VAR), data = .x),
        .options = furrr::furrr_options(
          globals = c("survfit", "update"),
          scheduling = 2
        ),
        .progress = TRUE
      ),
      # summarize K-M estimates at target event times
      km_summary = purrr::map(
        event_surv,
        ~ .get_surv_strata(.x, target_times)
      )
    ) |>
    tidyr::unnest(km_summary) |>
    mutate(VAR = as.integer(VAR))
  stopifnot(p_event |>
    dplyr::add_count(state, time, VAR) |>
    dplyr::filter(n > 1) |>
    nrow() == 0)

  # now, project K-M estimates onto test pop, for each subject in surv_vars group
  p_event_per_subject <- test_observed_events |>
    distinct(subject_id, !!!surv_vars) |>
    left_join(p_event_groups, by = dplyr::join_by(!!!surv_vars)) |>
    left_join(p_event, by = c("VAR"), multiple = "all", relationship = "many-to-many") |>
    dplyr::transmute(subject_id, state, time, p_event_free = surv) |>
    dplyr::arrange(state, time, subject_id)

  dplyr::lst(
    p_event, p_event_per_subject,
    prepared_observed_events, observed_events
  )
}


.summarize_ppsurv <- function(df, target_times, form = ~subject_id, p = NULL) {
  a <- survival::survfit(update(survival::Surv(time, is_event) ~ 1, form), data = df) |>
    summary(target_times)
  if ("strata" %in% names(a)) {
    res <- dplyr::as_tibble(dplyr::lst(
      time = a$time, surv = a$surv,
      strata = a$strata, cumhaz = a$cumhaz
    )) |>
      rgeco:::.tidy_km_strata()
  } else {
    if (is.null(a$surv)) {
      a <- survival::survfit(update(survival::Surv(time, is_event) ~ 1, form), data = df) |>
        summary(max(df$time))
      res <- dplyr::as_tibble(dplyr::lst(time = a$time, surv = a$surv, cumhaz = a$cumhaz)) |>
        dplyr::select(-time) |>
        tidyr::expand_grid(time = target_times)
    } else {
      res <- dplyr::as_tibble(dplyr::lst(time = a$time, surv = a$surv, cumhaz = a$cumhaz))
    }
  }
  if (!all(target_times %in% res$time)) {
    res <- res |>
      dplyr::bind_rows(res |>
        dplyr::select(-time, -surv, -cumhaz) |>
        dplyr::distinct() |>
        dplyr::mutate(surv = 0, cumhaz = max(res$cumhaz)) |>
        tidyr::expand_grid(time = target_times[!target_times %in% res$time]))
  }
  if (!is.null(p)) {
    p()
  }
  res
}

.summarize_ppsurv_old <- function(df, target_times, form = ~subject_id, p = NULL) {
  cov_syms <- .form_rhs_syms(form)
  df_groups <- df |>
    dplyr::mutate(all = 1) |>
    dplyr::group_by(all, !!!cov_syms) |>
    dplyr::summarize(n_total = dplyr::n())
  df_cumevents <- df |>
    dplyr::filter(is_event == 1) |>
    dplyr::mutate(all = 1) |>
    dplyr::group_by(all, !!!cov_syms, time) |>
    dplyr::summarize(cumevents = dplyr::n())
  if (nrow(df_cumevents) > 0) {
    df_cumevents <- df_cumevents |>
      dplyr::group_by(all, !!!cov_syms) |>
      dplyr::arrange(all, !!!cov_syms, time) |>
      dplyr::mutate(cumevents = purrr::accumulate(cumevents, .f = sum)) |>
      dplyr::select(all, !!!cov_syms, time, cumevents)
  }

  df_cumcensored <- df |>
    dplyr::filter(is_event == 0) |>
    dplyr::mutate(time = time + 0.1) |>
    dplyr::mutate(all = 1) |>
    dplyr::group_by(all, !!!cov_syms) |>
    dplyr::group_by(all, !!!cov_syms, time) |>
    dplyr::summarize(cumcensored = dplyr::n())
  if (nrow(df_cumcensored) > 0) {
    df_cumcensored <- df_cumcensored |>
      dplyr::group_by(all, !!!cov_syms) |>
      dplyr::arrange(all, !!!cov_syms, time) |>
      dplyr::mutate(cumcensored = purrr::accumulate(cumcensored, .f = sum)) |>
      dplyr::ungroup() |>
      dplyr::select(all, !!!cov_syms, time, cumcensored)
  }

  res <- df_groups |>
    tidyr::expand_grid(time = target_times) |>
    dplyr::bind_rows(df_cumcensored) |>
    dplyr::bind_rows(df_cumevents) |>
    dplyr::group_by(all, !!!cov_syms) |>
    dplyr::arrange(all, !!!cov_syms, time) |>
    tidyr::fill(cumevents, cumcensored, .direction = c("down")) |>
    dplyr::ungroup() |>
    dplyr::filter(time %in% target_times) |>
    dplyr::mutate_at(
      .vars = dplyr::vars(cumevents, cumcensored),
      .funs = ~ dplyr::if_else(is.na(.x), 0, .x)
    ) |>
    dplyr::mutate(
      n_at_risk = n_total - cumcensored,
      p_event = cumevents / n_at_risk,
      surv = 1 - p_event
    )
  if (!is.null(p)) {
    p()
  }
  res
}

.convert_str_to_form <- function(str) {
  if (!is.null(str) && length(str) > 0) {
    by_formula <- as.formula(str_c("~ ", str_c(str, collapse = "+")))
  } else {
    by_formula <- ~1
  }
  by_formula
}

#' Summarize paths using 'ppsurv'
#'
#' @export
#' @param pd A \code{\link{PathData}} object
#' @param target_times Times where to summarize
#' @param by Factor by which to summarize
#' @param split_by Factor by which to split
#' @param truncate_at_terminal_events Truncate paths at terminal events?
#' TODO: explain this better
#' @param use_old Use old version?
summarize_ppsurv <- function(pd, target_times, by = "subject_id",
                             split_by = by,
                             truncate_at_terminal_events = TRUE,
                             use_old = FALSE) {
  if (is.null(split_by)) {
    split_by <- c()
  }
  if (isTRUE(use_old)) {
    .surv_fun <- .summarize_ppsurv_old
  } else {
    .surv_fun <- .summarize_ppsurv
  }
  message("preparing data")
  # construct list of data frames for processing, by state w/wout split_by
  event_states <- pd$get_event_states()
  split_syms <- rlang::syms(split_by)
  by_formula <- .convert_str_to_form(setdiff(by, split_by))
  by_syms <- rlang::syms(setdiff(by, split_by))

  res_data <- pd$as_data_frame(
    truncate = truncate_at_terminal_events,
    covariates = unique(c(split_by, by))
  ) |>
    as_time_to_first_event(
      states = event_states,
      by = unique(c(split_by, by))
    ) |>
    mutate(grp = dense_rank(str_c(!!!split_syms, state))) |>
    dplyr::select(time, is_event, state, grp, !!!split_syms, !!!by_syms) |>
    dplyr::group_by(state, grp, !!!split_syms) |>
    tidyr::nest(data = c(time, is_event, !!!by_syms)) |>
    dplyr::ungroup()

  res_data$data |>
    rlang::set_names(res_data$grp) |>
    furrr::future_map(
      .surv_fun,
      target_times = target_times,
      form = by_formula,
      .progress = TRUE,
      .options =
        furrr::furrr_options(
          scheduling = 2,
          globals = c(
            ".summarize_ppsurv",
            ".summarize_ppsurv_old",
            ".surv_fun",
            ".form_rhs_syms"
          )
        )
    ) |>
    dplyr::bind_rows(.id = "grp") |>
    mutate(grp = as.integer(grp)) |>
    left_join(res_data |> distinct(grp, state, !!!split_syms))
}

.Cindex2 <- function(object, predicted, t_star = -1) {
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- predicted
    t_star0 <- t_star
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] ==
      1])))
    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    }
    vec_coxph <- predictSurvProb(obj, test_data, t_star0)
    object_coxph <- survival::Surv(test_data$time, test_data$status)
    object <- object_coxph
    predicted <- vec_coxph
  }
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- predicted
    t_star0 <- t_star
    distime <- obj$time.interest
    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    }
    med_index <- order(abs(distime - t_star0))[1]
    mat_rsf <- predict(obj, test_data)$survival
    vec_rsf <- mat_rsf[, med_index]
    object_rsf <- survival::Surv(test_data$time, test_data$status)
    object <- object_rsf
    predicted <- vec_rsf
  }
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- predicted
    t_star0 <- t_star
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] ==
      1])))
    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    }
    predicted <- predictSurvProb2survreg(
      obj, test_data,
      t_star0
    )
    object <- survival::Surv(test_data$time, test_data$status)
  }
  time <- object[, 1]
  status <- object[, 2]
  if (length(time) != length(status)) {
    stop("The lengths of time and status are not equal")
  }
  if (length(time) != length(predicted)) {
    stop("The lengths of time and predicted are not equal")
  }
  if (any(is.na(time) | is.na(status) | is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  permissible <- 0
  concord <- 0
  par_concord <- 0
  n <- length(time)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((time[i] < time[j] & status[i] == 0)) {
        next
      }
      if ((time[j] < time[i] & status[j] == 0)) {
        next
      }
      if (time[i] == time[j] & status[i] == 0 & status[j] ==
        0) {
        next
      }
      permissible <- permissible + 1
      if (time[i] != time[j]) {
        if ((time[i] < time[j] & predicted[i] < predicted[j])) {
          concord <- concord + 1
        } else if ((time[j] < time[i] & predicted[j] < predicted[i])) {
          concord <- concord + 1
        } else if (predicted[i] == predicted[j]) {
          par_concord <- par_concord + 0.5
        }
      }
      if (time[i] == time[j] & status[i] == 1 & status[j] ==
        1) {
        if (predicted[i] == predicted[j]) {
          concord <- concord + 1
        } else {
          par_concord <- par_concord + 0.5
        }
      }
      if (time[i] == time[j] & ((status[i] == 1 & status[j] ==
        0) | (status[i] == 0 & status[j] == 1))) {
        if ((status[i] == 1 & predicted[i] < predicted[j]) |
          (status[j] == 1 & predicted[j] < predicted[i])) {
          concord <- concord + 1
        } else {
          par_concord <- par_concord + 0.5
        }
      }
    }
  }
  C_index <- (concord + par_concord) / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}

#' More efficient version of SurvMetrics::Brier
.brier_score2 <- function(object, pre_sp, t_star, gt_times = NULL,
                          gt_tstar = NULL) {
  if (is.na(t_star)) {
    stop("Cannot calculate Brier Score at NA")
  }
  if (t_star <= 0) {
    stop("t_star must be strictly positive")
    # possible default: use median observed time
    # t_star <- median(object[, 1][object[, 2] == 1])
  }
  if (length(t_star) != 1) {
    stop("Brier Score can only be calculated at a single time point")
  }
  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }
  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(pre_sp))) {
    stop("The input probability vector cannot have NA")
  }
  if (length(object) != length(pre_sp)) {
    stop("The prediction survival probability and the survival object have different lengths")
  }
  time <- object[, 1]
  status <- object[, 2]
  t_order <- order(time)
  time <- time[t_order]
  status <- status[t_order]
  pre_sp <- pre_sp[t_order]
  if (is.null(gt_times)) {
    gt_times <- .Gt(object, time)
  } else {
    gt_times <- gt_times[t_order]
  }
  if (is.null(gt_tstar)) {
    gt_tstar <- .Gt(object, t_star)
  }
  sum_before_t <- 0
  sum_after_t <- 0
  # TODO: can be optimized because time is ordered
  index_events_before_tstar <- which(time < t_star & status == 1)
  index_censored_after_tstar <- which(time >= t_star)
  for (i in index_events_before_tstar) {
    sum_before_t <- sum_before_t + 1 / gt_times[i] * (pre_sp[i])^2
  }
  for (i in index_censored_after_tstar) {
    sum_after_t <- sum_after_t + 1 / gt_tstar * (1 - pre_sp[i])^2
  }
  BSvalue <- (sum_before_t + sum_after_t) / length(time)
  names(BSvalue) <- "Brier Score"
  return(round(BSvalue, 6))
}

.Gt <- function(object, timepoints) {
  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }
  if (missing(object)) {
    stop("The survival object is missing")
  }
  if (missing(timepoints)) {
    stop("The time is missing with no default")
  }
  if (any(timepoints < 0)) {
    stop("The timepoints must be positive")
  }
  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(timepoints))) {
    stop("Cannot calculate Gt at NA")
  }
  time <- object[, 1]
  status <- object[, 2]
  status0 <- ifelse(status == min(status), 1, 0)
  fit <- survival::survfit(survival::Surv(time, status0) ~ 1)
  res.sum <- survminer::surv_summary(fit)
  res.sum <- na.omit(res.sum) |>
    dplyr::bind_rows(tidyr::expand_grid(time = timepoints) |>
      dplyr::mutate(index = seq_along(timepoints))) |>
    dplyr::arrange(time, surv) |>
    # use stepwise censoring function instead of interpolating as done in the
    # origin `Gt` function
    tidyr::fill(surv, .direction = "downup") |>
    dplyr::mutate(surv = dplyr::if_else(is.na(surv) & time == min(time), 1, surv))
  if (all(timepoints %in% res.sum$time)) {
    Gvalue <- res.sum |>
      dplyr::filter(!is.na(index)) |>
      dplyr::arrange(index) |>
      dplyr::pull(surv)
  } else {
    stop("Error computing censoring function: not all times represented")
  }
  if (any(is.na(Gvalue) | Gvalue <= 0)) {
    stop("Error computing censoring function: some values are invalid")
  }
  return(Gvalue)
}

.integrate_brier_score <- function(bs_estimate, by = c()) {
  by <- rlang::syms(by)
  bs_estimate |>
    dplyr::group_by(state, method, Event, !!!by) |>
    dplyr::arrange(time) |>
    mutate(
      interval = dplyr::lead(time, n = 1, order_by = time) - time,
      ibs = purrr::accumulate(brier_score * interval, sum)
    ) |>
    dplyr::ungroup()
}

.common_subjects <- function(pd, pred_pd = NULL,
                             ppsurv = NULL, sub_ids_char,
                             sub_ids_char_training = NULL,
                             .label = "result") {
  # construct list of subject ids in common across all inputs
  obs_subs <- unique(pd$subject_df$subject_id)

  if (!is.null(pred_pd)) {
    pred_subs <- unique(pred_pd$subject_df$subject_id)
  } else if (!is.null(ppsurv)) {
    pred_subs <- unique(ppsurv$subject_id)
  } else {
    pred_subs <- obs_subs
  }
  if (is.null(sub_ids_char)) {
    sub_ids_char <- intersect(obs_subs, pred_subs)
    if (length(obs_subs) != length(sub_ids_char)) {
      message(glue::glue("Note: computing {str_c(.label, collapse = '+')} for the {length(sub_ids_char)} subjects in common between observed (n={length(obs_subs)}) and predicted (n={length(pred_subs)}) data."))
    }
  } else {
    sub_ids <- intersect(intersect(obs_subs, pred_subs), sub_ids_char)
    if (length(sub_ids_char) != length(sub_ids)) {
      message(glue::glue("Note: computing {str_c(.label, collapse = '+')} for the {length(sub_ids)} subjects in common between observed (n={length(obs_subs)}), predicted (n={length(pred_subs)}), and user-provided (n={length(sub_ids_char)}) data."))
    }
    sub_ids_char <- sub_ids
  }
  sub_ids_char
}

.common_subjects_training <- function(pd, sub_ids_char,
                                      sub_ids_char_training = NULL,
                                      .label = "K-M expectation") {
  # construct list of subject ids in common across all inputs
  obs_subs <- unique(pd$subject_df$subject_id)
  if (is.null(sub_ids_char_training)) {
    message("Note: training sub_ids not provided. Assuming training and test sub_ids are the same.")
    sub_ids_char_training <- sub_ids_char
  } else {
    sub_ids <- intersect(obs_subs, sub_ids_char_training)
    if (length(sub_ids_char) != length(sub_ids)) {
      message(glue::glue("Note: Using the {length(sub_ids)} subjects in common between observed (n={length(obs_subs)}) and user-provided (n={length(sub_ids_char_training)}) training IDs to estimate {str_c(.label, collapse = '+')} in test set (n={length(sub_ids_char)})."))
    }
    sub_ids_char_training <- sub_ids
  }
  sub_ids_char_training
}

compute_scores <- function(pd, pred_pd = NULL,
                           target_times = NULL,
                           sub_ids_char = NULL,
                           sub_ids_char_training = NULL,
                           km_formula = ~1,
                           ppsurv = NULL,
                           keep_gt_info = FALSE,
                           which = c("brier_score", "cindex"),
                           by = c()) {
  by_syms <- rlang::syms(by)
  if (is.null(target_times) && !is.null(ppsurv)) {
    target_times <- unique(ppsurv$time)
  } else if (is.null(target_times)) {
    stop("target_times are required unless ppsurv is pre-computed.")
  }
  target_times <- target_times[target_times > 0]
  sub_ids_char <- .common_subjects(pd, pred_pd, ppsurv, sub_ids_char,
    .label = str_c(which, sep = " and ")
  )
  if (!is.null(sub_ids_char_training)) {
    sub_ids_char_training <- .common_subjects_training(
      pd,
      sub_ids_char = sub_ids_char,
      sub_ids_char_training = sub_ids_char_training,
      .label = str_c(which, sep = " and ")
    )
  }
  futile.logger::flog.info("Starting preparation of observed data ...")
  observed <- .prepare_observed_data(pd, target_times,
    sub_ids_char = sub_ids_char,
    sub_ids_char_training = sub_ids_char_training,
    formula = km_formula
  )
  futile.logger::flog.info("Finished preparation of observed data.")
  if (is.null(ppsurv) && !is.null(pred_pd)) {
    predicted <- pred_pd$filter(subject_ids_keep = sub_ids_char) |>
      summarize_ppsurv(
        target_times = target_times,
        by = unique(c("subject_id", by))
      ) |>
      mutate(p_event_free = surv)
  } else if (!is.null(ppsurv)) {
    predicted <- ppsurv |>
      dplyr::filter(
        subject_id %in% sub_ids_char,
        time %in% target_times
      ) |>
      mutate(p_event_free = surv)
  } else {
    predicted <- NULL
  }
  futile.logger::flog.info("Finished preparation of predicted data ...")

  # observed data for each target time & state
  observed_data_per_state <- observed$prepared_observed_events |>
    dplyr::select(-event_data)

  # expected survival probability according to K-M estimate
  km_prediction <- observed$p_event_per_subject |>
    dplyr::group_by(state, time) |>
    tidyr::nest(predicted = c(subject_id, p_event_free)) |>
    dplyr::ungroup() |>
    mutate(
      p_event_free = map(predicted, ~ dplyr::arrange(.x, subject_id) |>
        dplyr::pull(p_event_free)),
    ) |>
    dplyr::select(-predicted)

  # fix a rare but important edge case, if each lengths aren't the same
  equal_lengths <- (km_prediction |>
    mutate(n_p = map_int(p_event_free, length)) |>
    distinct(n_p) |> nrow() == 1)
  if (isFALSE(equal_lengths)) {
    km_prediction <- observed$p_event_per_subject |>
      tidyr::expand(subject_id, time, state) |>
      left_join(
        observed$p_event_per_subject |>
          dplyr::select(time, state, subject_id, p_event_free),
        by = dplyr::join_by(time, state, subject_id)
      ) |>
      mutate(p_event_free = if_else(is.na(p_event_free), 0, p_event_free)) |>
      dplyr::group_by(state, time) |>
      tidyr::nest(predicted = c(subject_id, p_event_free)) |>
      dplyr::ungroup() |>
      mutate(
        p_event_free = map(
          predicted,
          ~ dplyr::arrange(.x, subject_id) |> dplyr::pull(p_event_free)
        ),
      ) |>
      dplyr::select(-predicted)
  }

  km_scores <- NULL
  if ("brier_score" %in% which) {
    futile.logger::flog.info("Starting computation of brier score on K-M estimates ...")
    km_scores <-
      km_prediction |>
      dplyr::inner_join(observed_data_per_state, by = dplyr::join_by(state)) |>
      dplyr::mutate(
        score = furrr::future_pmap_dbl(
          dplyr::lst(
            object = st_surv,
            pre_sp = p_event_free,
            gt_times,
            t_star = time
          ), .brier_score2,
          .progress = TRUE,
          .options = furrr::furrr_options(
            scheduling = 2,
            globals = c(".brier_score2", ".Gt", "survfit")
          )
        ),
        score_name = "brier_score",
        method = "km"
      ) |>
      dplyr::bind_rows(km_scores)
  }
  if ("cindex" %in% which) {
    futile.logger::flog.info("Starting computation of cindex on K-M estimates ...")
    km_scores <-
      km_prediction |>
      dplyr::inner_join(observed_data_per_state, by = dplyr::join_by(state)) |>
      dplyr::mutate(
        score = furrr::future_pmap_dbl(
          dplyr::lst(
            object = st_surv,
            predicted = p_event_free,
            t_star = time
          ), purrr::possibly(.Cindex2, otherwise = NA_real_),
          .progress = TRUE,
          .options = furrr::furrr_options(scheduling = 2, globals = c(".Cindex2"))
        ),
        score_name = "cindex",
        method = "km"
      ) |>
      dplyr::bind_rows(km_scores)
  }

  msm_scores <- NULL
  if (!is.null(predicted)) {
    # scores for multistate model (predictions)
    futile.logger::flog.info("Starting preparation of msm data ...")
    msm_prediction <- predicted |>
      tidyr::expand(subject_id, time, state, !!!by_syms) |>
      left_join(predicted |> dplyr::select(time, state, subject_id, p_event_free, !!!by_syms),
        by = dplyr::join_by(time, state, subject_id, !!!by_syms)
      ) |>
      mutate(p_event_free = if_else(is.na(p_event_free), 0, p_event_free)) |>
      tidyr::nest(predicted = c(subject_id, p_event_free)) |>
      mutate(
        predicted = map(predicted, dplyr::arrange, subject_id),
        p_event_free = map(predicted, dplyr::pull, p_event_free)
      ) |>
      dplyr::select(-predicted)

    if ("brier_score" %in% which) {
      futile.logger::flog.info("Starting computation of brier score for msm data ...")
      msm_scores <- msm_prediction |>
        inner_join(observed_data_per_state, by = dplyr::join_by(state)) |>
        dplyr::mutate(
          score = furrr::future_pmap_dbl(
            dplyr::lst(
              object = st_surv,
              pre_sp = p_event_free,
              gt_times, t_star = time
            ),
            purrr::possibly(.brier_score2, otherwise = NA_real_),
            .progress = TRUE,
            .options = furrr::furrr_options(scheduling = 2, globals = c(".brier_score2", ".Gt", "survfit"))
          ),
          score_name = "brier_score",
          method = "msm"
        ) |>
        dplyr::bind_rows(msm_scores)
    }
    if ("cindex" %in% which) {
      futile.logger::flog.info("Starting computation of cindex for msm data ...")
      msm_scores <- msm_prediction |>
        inner_join(observed_data_per_state, by = dplyr::join_by(state)) |>
        dplyr::mutate(
          score = furrr::future_pmap_dbl(
            dplyr::lst(
              object = st_surv,
              predicted = p_event_free,
              t_star = time
            ),
            purrr::possibly(.Cindex2, otherwise = NA_real_),
            .progress = TRUE,
            .options = furrr::furrr_options(scheduling = 2, globals = c(".Cindex2"))
          ),
          score_name = "cindex",
          method = "msm"
        ) |>
        dplyr::bind_rows(msm_scores)
    }
  }

  futile.logger::flog.info("Starting formatting for output ...")
  bs_estimates <- dplyr::lst(
    km = km_scores,
    msm = msm_scores
  ) |>
    dplyr::bind_rows() |>
    tidyr::spread(score_name, score) |>
    dplyr::mutate(Event = pd$state_names[state])

  if (isFALSE(keep_gt_info)) {
    bs_estimates <- bs_estimates |>
      dplyr::select(-p_event_free, -st_surv, -gt_times)
  }
  if ("brier_score" %in% which) {
    bs_estimates <- bs_estimates |>
      .integrate_brier_score(by = by)
  }
  futile.logger::flog.info("Completed compute_scores.")
  bs_estimates
}


compute_d_calibration <- function(pd, pred_pd,
                                  sub_ids_char = NULL,
                                  truncate_at_terminal_events = FALSE) {
  # construct list of subject ids in common across all inputs
  sub_ids_char <- .common_subjects(pd, pred_pd, sub_ids_char,
    .label = "D-calibration"
  )

  observed_events <- pd$as_time_to_first_event(
    truncate = truncate_at_terminal_events
  ) |>
    dplyr::filter(subject_id %in% sub_ids_char)
  observed_times <- unique(observed_events$time) |> sort()
  if (!is.null(pred_pd)) {
    predicted <- .prepare_predicted_data(pred_pd,
      target_times = observed_times,
      sub_ids_char = sub_ids_char, covs = c("subject_id")
    )
  }

  # get Pr(survival) for each observed event time and each state
  merged_data <- observed_events |>
    dplyr::select(subject_id, state, time, is_event, Event) |>
    left_join(
      predicted$p_event_per_subject |>
        dplyr::select(subject_id, state, time, proportion, p_event_free),
      by = c("subject_id", "state", "time")
    )

  all_qtiles <- tibble(
    p_qtile_lb = seq(from = 0, to = 0.9, by = 0.1),
    p_qtile_ub = seq(from = 0.1, to = 1, by = 0.1),
    value = seq(from = 0.05, to = 0.95, by = 0.1)
  ) |>
    mutate(p_qtile = cut_interval(value, length = 0.1, ordered_result = TRUE)) |>
    dplyr::select(-value)

  qtiles_observed <- merged_data |>
    dplyr::filter(is_event == 1) |>
    cross_join(all_qtiles) |>
    dplyr::filter(p_event_free <= p_qtile_ub, p_event_free > p_qtile_lb) |>
    dplyr::group_by(p_qtile, Event, state) |>
    dplyr::tally() |>
    dplyr::ungroup() |>
    dplyr::arrange(Event, state, p_qtile)

  qtiles_censored <- merged_data |>
    dplyr::filter(is_event == 0) |>
    cross_join(all_qtiles) |>
    dplyr::filter(p_event_free >= p_qtile_lb | p_event_free > p_qtile_ub) |>
    mutate(n = case_when(
      p_event_free >= p_qtile_lb ~ 1 - (p_qtile_lb / p_event_free),
      TRUE ~ 1 / (10 * p_event_free)
    )) |>
    dplyr::group_by(Event, state, p_qtile) |>
    summarize(n = sum(n)) |>
    dplyr::ungroup() |>
    dplyr::arrange(Event, state, p_qtile)

  d_calibration <- list(censored = qtiles_censored, uncensored = qtiles_observed) |>
    dplyr::bind_rows(.id = "type")
}

compute_one_calibration <- function(pd, pred_pd = NULL, ppsurv = NULL,
                                    target_times = seq(from = 0, to = 3 * 365.25, by = 365.25),
                                    sub_ids_char = NULL,
                                    sub_ids_char_training = NULL,
                                    truncate_at_terminal_events = FALSE) {
  # construct list of subject ids in common across all inputs
  sub_ids_char <- .common_subjects(
    pd = pd, pred_pd = pred_pd,
    ppsurv = ppsurv,
    sub_ids_char = sub_ids_char,
    .label = "1-calibration"
  )
  if (!is.null(sub_ids_char_training)) {
    sub_ids_char <- .common_subjects_training(pd,
      sub_ids_char = sub_ids_char,
      sub_ids_char_training = sub_ids_char_training,
      .label = "1-calibration"
    )
  }
  observed_events <- pd$filter(
    subject_ids_keep = sub_ids_char
  )$as_data_frame(
    covariates = c("subject_id")
  ) |>
    dplyr::filter(is_event == 1) |>
    dplyr::group_by(subject_id, state) |>
    dplyr::filter(time == min(time)) |>
    transmute(subject_id, state, observed_event_time = time)
  state_decode <- tibble(
    state = pd$get_event_states(),
    Event = pd$get_event_state_names()
  )
  if (is.null(ppsurv) && !is.null(pred_pd)) {
    predicted <- pred_pd$filter(subject_ids_keep = sub_ids_char) |>
      summarize_ppsurv(
        target_times = target_times,
        by = unique(c("subject_id", by)),
        split_by = "subject_id",
        truncate_at_terminal_events = truncate_at_terminal_events
      ) |>
      mutate(
        p_event_free = surv,
        proportion = 1 - surv
      )
  } else if (!is.null(ppsurv)) {
    predicted <- ppsurv |>
      dplyr::filter(
        subject_id %in% sub_ids_char,
        time %in% target_times
      ) |>
      mutate(
        p_event_free = surv,
        proportion = 1 - surv
      )
  } else {
    warning("Predicted data not provided! Returning NULL")
    return(NULL)
  }

  # get Pr(survival) for each state and target_time
  merged_data <- predicted |>
    dplyr::ungroup() |>
    mutate(n_subjects = n_distinct(subject_id)) |>
    dplyr::filter(time %in% target_times, time > 0, state > 1) |>
    left_join(observed_events,
      by = c("subject_id", "state"),
      multiple = "all", relationship = "many-to-one"
    ) |>
    mutate(event_observed = case_when(
      is.na(observed_event_time) ~ 0,
      observed_event_time > time ~ 0,
      observed_event_time <= time ~ 1
    )) |>
    dplyr::group_by(time, state) |>
    mutate(psurv_bin = ntile(proportion, n = 10)) |>
    dplyr::group_by(time, state, psurv_bin) |>
    summarize(
      predicted_events = sum(proportion),
      observed_events = sum(event_observed),
      n = n(),
      p = mean(proportion)
    ) |>
    mutate(test_stat = (
      (observed_events - predicted_events)^2) / (predicted_events * (1 - p))) |>
    tidyr::gather(type, value, observed_events, predicted_events) |>
    left_join(state_decode)

  test_stat <- merged_data |>
    distinct(time, state, psurv_bin, test_stat) |>
    dplyr::group_by(time, state) |>
    summarise(
      HLstat = sum(test_stat),
      p_value = pchisq(HLstat, 10 - 2, lower.tail = FALSE)
    )
  merged_data |>
    dplyr::select(-test_stat, -p, -n) |>
    left_join(test_stat)
}

write_score_to_file <- function(bs_estimates, params, cache_name,
                                dir = here::here("brier_scores")) {
  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir, recurse = T)
  }
  file <- fs::path(dir, fs::path_ext_set(cache_name, ".csv"))
  if (fs::file_exists(file)) {
    append <- TRUE
    if (interactive()) {
      if (askYesNo("File already exists. Replace?")) {
        fs::file_delete(file)
        append <- FALSE
      }
    }
  } else {
    append <- FALSE
  }

  # prepare run info as a row in a dataframe
  run_vars <- vars(
    project, project_version, sample_n, events,
    iter_warmup, iter_sampling,
    chains, NK, cache_name
  )
  params$covariates <- NULL
  run_data <- tibble(!!!params) |>
    mutate(cache_name = cache_name) |>
    dplyr::select(!!!run_vars)

  # prepare IBS summary per event type
  result_vars <- vars(n_subjects, time, Event, interval, km, msm)
  result_summary <- bs_estimates |>
    dplyr::select(Event, method, time, ibs, n_subjects, interval) |>
    drop_na() |>
    dplyr::filter(time %in% c(max(time), seq(from = 0, to = 3 * 365.25, by = 365.25))) |>
    spread(method, ibs) |>
    dplyr::select(!!!result_vars)

  results <- result_summary |>
    bind_cols(run_data) |>
    mutate(date = lubridate::today()) |>
    dplyr::select(date, !!!run_vars, !!!result_vars)

  write_csv(results, file = file, append = append, quote = "all")
}


# Plot
plot_bs <- function(bs_est) {
  bs_est |>
    ggplot(aes(x = time, y = brier_score, color = method)) +
    facet_wrap(~Event, scale = "free_y") +
    geom_line() +
    scale_y_continuous("Brier Score (MSE)") +
    theme(legend.position = "top")
}

# Is a pair comparable?
is_comparable <- function(y1, y2, censored) {
  isFALSE(censored) && (y1 < y2)
}

# Concordance index for right censored data
concordance_index <- function(t_true, t_pred, censored) {
  checkmate::assert_numeric(t_true, min.len = 1)
  L <- length(t_true)
  checkmate::assert_numeric(t_pred, len = L)
  checkmate::assert_logical(censored, len = L)
  st <- survival::Surv(t_true, !censored)
  SurvMetrics::Cindex(st, t_pred)
}

# Concordance index given path data
c_index <- function(pd_obs, pd_pred, state_idx, kmfit, subject_ids, t_eval) {
  message("Extracting (state ", state_idx, ")")
  fet <- first_events(pd_obs, pd_pred, state_idx, kmfit, subject_ids, t_eval)
  df_et <- fet$df
  sn <- pd_obs$state_names[state_idx]
  fet$km_plot <- fet$km_plot + ggtitle(sn)
  message("Computing concordance index")
  out <- list(
    ci_t_msm = concordance_index(df_et$time, df_et$t_msm, !df_et$observed),
    ci_p_msm = concordance_index(df_et$time, 1 - df_et$p_msm, !df_et$observed),
    ci_p_km = concordance_index(df_et$time, 1 - df_et$p_km, !df_et$observed),
    state_name = sn
  )

  # Table format
  table <- data.frame(
    method = c("km", "msm"),
    cindex = c(out$ci_p_km, out$ci_p_msm)
  )
  table$Event <- sn

  # Return
  out$table <- table
  c(fet, out)
}

# Kaplan-Meier fit
km_fit <- function(pd, state_idx, subject_ids) {
  checkmate::assert_number(state_idx)
  covs <- c("subject_id", "first_dose_amount", "dose_adjustment")
  df_obs <- pd$as_data_frame(covariates = covs) |>
    dplyr::filter(subject_id %in% subject_ids)
  te_obs <- time_to_event(df_obs, state_idx) |> dplyr::select(-path_id)
  te_obs <- add_dose_arm(te_obs, df_obs)
  survival::survfit(survival::Surv(time, observed) ~ dose_arm, te_obs)
}

# Predict survival probability for each dose arm using previous km fit
km_predict <- function(kmfit, time) {
  checkmate::assert_number(time)
  a <- summary(kmfit, times = time)
  da <- sapply(
    strsplit(as.vector(a$strata), split = "=", fixed = TRUE),
    function(x) x[2]
  )
  data.frame(dose_arm = da, p_km = 1 - a$surv)
}

# Add dose arm to df
add_dose_arm <- function(te_obs, df) {
  df$dose_arm <- paste0(df$dose_adjustment, "-", df$first_dose_amount)
  df <- df |>
    dplyr::select(subject_id, dose_arm) |>
    dplyr::ungroup() |>
    dplyr::group_by(subject_id) |>
    dplyr::slice(1)
  te_obs |> left_join(df, by = "subject_id")
}

# First event times and event probabilities given path data
first_events <- function(pd_obs, pd_pred, state_idx, kmfit, subject_ids, tev) {
  # Filter subjects
  covs <- c("subject_id", "dose_adjustment", "first_dose_amount")
  df_obs <- pd_obs$as_data_frame(covariates = covs) |>
    dplyr::filter(subject_id %in% subject_ids)
  df_pred <- pd_pred$as_data_frame(covariates = covs) |>
    dplyr::filter(subject_id %in% subject_ids)

  # Times to first event
  te_pred <- time_to_event(df_pred, state_idx) |> dplyr::select(-path_id)
  te_med <- te_pred |>
    dplyr::group_by(subject_id) |>
    summarize(t_msm = mean(time))
  te_prop <- te_pred |>
    dplyr::group_by(subject_id) |>
    summarize(p_msm = mean(observed))
  te_obs <- time_to_event(df_obs, state_idx) |> dplyr::select(-path_id)
  te_obs <- add_dose_arm(te_obs, df_obs)

  # Kaplan-Maier
  df_km <- km_predict(kmfit, tev)
  kms <- survminer::surv_summary(kmfit, te_obs)
  sp <- ggplot(kms, aes(x = time, y = surv, color = dose_arm)) +
    geom_line()

  # Return
  list(
    df = te_obs |>
      left_join(te_med, by = "subject_id") |>
      left_join(te_prop, by = "subject_id") |>
      left_join(df_km, by = "dose_arm"),
    km_plot = sp,
    km_fit = km_fit,
    df_km = df_km
  )
}



# Plot concordance index
plot_cindex <- function(ci, name, nrow = NULL, ncol = NULL) {
  plots <- list()
  if (name == "ci_t_msm") {
    ylab <- "t_msm"
  } else if (name == "ci_p_msm") {
    ylab <- "p_msm"
  } else if (name == "ci_p_km") {
    ylab <- "p_km"
  } else {
    stop("invalid option")
  }
  for (j in 1:length(ci)) {
    res <- ci[[j]]
    rn <- res[[name]]
    df <- res$df |> dplyr::arrange(observed)
    st <- paste0("ci = ", round(rn, 4))
    plt <- df |> ggplot(aes(x = time, y = !!sym(ylab), color = observed)) +
      geom_point(alpha = 0.95, size = 0.4)
    plots[[j]] <- plt + ggtitle(names(ci)[j], subtitle = st) +
      theme(legend.position = "top") +
      scale_color_brewer(palette = 7, type = "seq", direction = 1) +
      xlab("") + ylab("") + theme(
        legend.position = "none", text = element_text(size = 8)
      )
  }
  plt <- ggdplyr::arrange(
    plotlist = plots, legend.grob = get_legend(plots[[1]]),
    nrow = nrow, ncol = ncol
  )
  annotate_figure(plt, top = paste0("Concordance index (", ylab, ")"))
}

# Plot multiple concordance indices
plot_ci_multi <- function(ci) {
  names(ci) <- sapply(ci, function(x) x$state_name)
  ci_mean <- sapply(ci, function(x) x$ci_mean)
  ci_prob <- sapply(ci, function(x) x$ci_prob)
  plt_1 <- plot_cindex(ci, name = "ci_p_km", nrow = 1)
  plt_2 <- plot_cindex(ci, name = "ci_p_msm", nrow = 1)
  ggdplyr::arrange(plt_1, plt_2, ncol = 1, nrow = 2)
}
