# Proportion of paths that has a given event occur at least once
p_event <- function(paths_df, state_names = NULL) {
  nonevent_states <- c("Censor", "Randomization")
  n_paths <- length(unique(paths_df$path_id))
  paths_df <- paths_df %>%
    arrange(path_id, time) %>%
    group_by(state) %>%
    filter(is_event == 1)
  df_empty <- data.frame(state_char = state_names, state = 1:length(state_names))
  df <- paths_df %>%
    summarise(p_paths = dplyr::n_distinct(path_id) / n_paths)
  df <- df %>% filter(state != 0)
  df <- df_empty %>% left_join(df, by = join_by(state))
  df$p_paths[is.na(df$p_paths)] <- 0
  df
  df %>% filter(!(state_char %in% nonevent_states))
}

# Above but for each subject separately
p_event_by_subject <- function(pd, draw_idx = NULL, rep_idx = NULL) {
  pd_df <- pd$as_data_frame(covariates = "subject_id")
  state_names <- pd$state_names
  if (!is.null(rep_idx)) {
    pd_df <- pd_df[which(pd_df$rep_idx == rep_idx), ]
  }
  if (!is.null(draw_idx)) {
    pd_df <- pd_df[which(pd_df$draw_idx == draw_idx), ]
  }
  a <- unique(pd$subject_df$subject_id)
  list_of_dfs <- list()
  J <- length(a)
  for (j in seq_len(J)) {
    df <- pd_df[which(pd_df$subject_id == a[j]), ]
    list_of_dfs[[j]] <- p_event(df, state_names)
  }
  names(list_of_dfs) <- a
  out <- stack_list_of_dfs(list_of_dfs)
  out$subject_id <- out$name
  out %>% select(-name)
}

# Stack list of data frames to one data frame
# Could probably be replaced by dplyr::bind_rows()
stack_list_of_dfs <- function(list_of_dfs) {
  bindfun <- function(x) {
    cbind(list_of_dfs[[x]], name = x)
  }
  do.call(rbind, lapply(names(list_of_dfs), bindfun))
}


# p_event by variable (list of pathdata)
.p_event_by_vars_list <- function(pd_list, vars = "first_dose_amount") {
  sym_vars <- rlang::syms(vars)
  pd_df <- pd_list |>
    map_dfr(~ .x$as_data_frame(covariates = vars), .id = "_sample") |>
    mutate(path_id = str_c(`_sample`, path_id)) |>
    mutate(`_val` = str_c(!!!sym_vars, sep = ":"))
  state_names <- pd_list[[1]]$state_names
  vals <- pd_df |>
    distinct(`_val`) |>
    pull(`_val`) |>
    sort()

  df <- vals |>
    set_names() |>
    purrr::map(~ dplyr::filter(pd_df, `_val` == .x)) |>
    furrr::future_map_dfr(~ p_event(.x, state_names), .id = "_val")

  df |>
    left_join(pd_df |> distinct(`_val`, !!!sym_vars)) |>
    select(-`_val`)
}

# p_event by variable
.p_event_by_vars <- function(pd, vars = "first_dose_amount",
                             draw_idx = NULL, rep_idx = NULL) {
  pd_df <- pd$filter(draw_ids_keep = draw_idx, rep_ids_keep = rep_idx)$as_data_frame(covariates = vars)
  sym_vars <- rlang::syms(vars)
  pd_df <- pd_df |>
    mutate(`_val` = str_c(!!!sym_vars, sep = ":"))
  vals <- pd_df |>
    distinct(`_val`) |>
    pull(`_val`) |>
    sort()
  state_names <- pd$state_names

  df <- vals |>
    set_names() |>
    purrr::map(~ dplyr::filter(pd_df, `_val` == .x)) |>
    furrr::future_map_dfr(~ p_event(.x, state_names),
      .id = "_val"
    )

  df |>
    left_join(pd_df |> distinct(`_val`, !!!sym_vars), by = join_by(`_val`)) |>
    select(-`_val`)
}

# Function for transforming format (counterfactual paths)
cf_p_event_by_subject <- function(list_of_pds, draw_idx = NULL, rep_idx = NULL) {
  ff <- function(x) {
    p_event_by_subject(x, draw_idx, rep_idx)
  }
  pes <- lapply(list_of_pds, ff)
  df_cfs <- stack_list_of_dfs(pes)
  f1 <- function(x) as.numeric(x[1])
  df_cfs$cf_dose <- sapply(strsplit(df_cfs$name, split = " ", fixed = T), f1)
  df_cfs
}

# Study variation over reps (counterfactual paths)
cf_p_event_by_subject_by_rep <- function(list_of_pds, n_repeats) {
  df <- NULL
  for (j in seq_len(n_repeats)) {
    df_j <- cf_p_event_by_subject(list_of_pds, rep_idx = j)
    df_j$rep_idx <- j
    df <- rbind(df, df_j)
  }
  df %>%
    group_by(state_char, cf_dose) %>%
    summarize(
      p_mean = mean(p_paths), p_std = sd(p_paths),
      .groups = "drop"
    )
}

# Study variation over posterior draws (counterfactual paths)
cf_p_event_by_subject_by_draw <- function(list_of_pds, n_draws) {
  df <- NULL
  for (j in seq_len(n_draws)) {
    df_j <- cf_p_event_by_subject(list_of_pds, draw_idx = j)
    df_j$draw_idx <- j
    df <- rbind(df, df_j)
  }
  df %>%
    group_by(state_char, cf_dose) %>%
    summarize(
      p_mean = mean(p_paths), p_std = sd(p_paths),
      .groups = "drop"
    )
}

# Plot result of cf_p_event_by_subject_by_rep for one subject
plot_cf_dose_effect <- function(p_cf, n_repeats) {
  n_draws <- length(unique((p_cf[[1]]$link_df$draw_idx)))
  message("Computing p_event by rep")
  df <- cf_p_event_by_subject_by_rep(p_cf, n_repeats)
  message("Computing p_event by draw")
  df2 <- cf_p_event_by_subject_by_draw(p_cf, n_draws)

  ggplot(df, aes(
    x = cf_dose, y = p_mean, ymin = p_mean - p_std, ymax = p_mean + p_std
  )) +
    geom_errorbar(width = 4) +
    geom_line() +
    geom_point() +
    geom_errorbar(data = df2, color = "red", width = 4) +
    geom_line(data = df2, color = "red", lty = 2) +
    facet_wrap(. ~ state_char, scales = "free_y", nrow = 1) +
    xlab("Dose (mg)") +
    ylab("3-year event rate")
}

# Proportion of paths that had a given event until given time in a grid
compute_proportion_up_to_time <- function(data, time_grid, event_states, byvars = vars()) {
  all_groups <- data |>
    distinct(!!!byvars) |>
    expand(!!!byvars) |>
    expand_grid(
      time = time_grid,
      state = event_states
    )
  if (nrow(data) == 0) {
    return(NULL)
  }
  results_pre <- data |>
    group_by(!!!byvars) |>
    mutate(total_n_paths = dplyr::n_distinct(path_id)) |>
    ungroup() |>
    filter(is_event == 1) |>
    group_by(!!!byvars, state, path_id) |>
    # keep first record where path observed in state
    filter(time == min(time)) |>
    group_by(!!!byvars, state) |>
    arrange(time) |>
    mutate(
      n_paths = 1,
      n_paths_visited = purrr::accumulate(n_paths, sum),
      p_paths = n_paths_visited / total_n_paths
    ) |>
    ungroup() |>
    distinct(!!!byvars, time, state, p_paths)
  results <- results_pre |>
    bind_rows(all_groups) |>
    group_by(!!!byvars, state) |>
    arrange(time) |>
    tidyr::fill(p_paths, .direction = "down") |>
    ungroup() |>
    mutate(proportion = if_else(is.na(p_paths), 0, p_paths)) |>
    distinct(!!!byvars, state, time, proportion) |>
    filter(time %in% time_grid)

  dplyr::as_tibble(results)
}

# Estimate probability of no event until time t
estimate_noevent_probability <- function(data, time_points) {
  p <- compute_proportion_up_to_time(data, time_points)
  p$p_noevent <- 1 - p$proportion
  p$proportion <- NULL
  p
}

# Time to first time state_idx is observed in observed paths
# Could be replaced by PathData$as_time_to_first_event()
time_to_event <- function(df, state_idx) {
  df %>%
    group_by(path_id, subject_id) %>%
    summarize(
      time = ifelse(any(state == state_idx),
        min(time[state == state_idx]), max(time)
      ),
      observed = any(state == state_idx),
      .groups = "drop"
    )
}
