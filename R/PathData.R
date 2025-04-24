# util
pathdata_size_per_row <- function(pd) {
  object.size(pd$df) / nrow(pd$df)
}

# util
check_columns <- function(df, needed_columns) {
  if (!(all(needed_columns %in% colnames(df)))) {
    message("found = {", paste(colnames(df), collapse = ", "), "}")
    message("required = {", paste(needed_columns, collapse = ", "), "}")
    stop("some needed columns are missing from df")
  }
}

# Initialize data from single observational path data frame that has
# One path per subject
create_pathdata <- function(df, covs, ...) {
  draw_idx <- 1
  rep_idx <- 1
  df$path_id <- paste0(df$subject_id, draw_idx, rep_idx)
  df$draw_idx <- 1
  df$rep_idx <- 1
  subject_df <- df |>
    group_by(subject_id) |>
    slice(1)
  link_df <- df |>
    group_by(path_id) |>
    slice(1)
  PathData$new(subject_df, df, link_df, covs = covs, ...)
}


#' Path data class (R6 class)
#'
#' @export
#' @field subject_df Data frame with one row per subject. Must have
#' \code{subject_id} and all covariates as columns.
#' @field path_df Data frame of actual paths. Must have \code{path_id},
#' \code{state}, \code{time}, and \code{is_event} as columns.
#' @field link_df Links the path and subject data frames. Must have
#' \code{path_id}, \code{draw_idx}, \code{rep_idx}, and \code{subject_id} as
#' columns.
#' @field covs Covariate column names.
#' @field state_names Names of the states.
#' @field state_types State types.
#' @field terminal_states Terminal states.
#' @field initial_states Initial states.
#' @field censor_states Censoring states.
#' @field dt TODO
PathData <- R6::R6Class(
  classname = "PathData",
  public = list(
    subject_df = NULL, # must have subject_id and all covariates
    path_df = NULL, # must have
    link_df = NULL, #
    covs = NULL,
    state_names = NULL,
    state_types = NULL,
    terminal_states = NULL,
    initial_states = NULL,
    censor_states = NULL,
    dt = NULL,

    #' Initialize
    #' @param subject_df Subjects data frame.
    #' @param path_df Paths data frame.
    #' @param link_df Link data frame.
    #' @param covs Covariate column names.
    #' @param check_order Check order of paths?
    #' @param state_names Names of the states.
    #' @param terminal_states Terminal states.
    #' @param initial_states Initial states.
    #' @param censor_states Censoring states.
    initialize = function(subject_df, path_df, link_df,
                          state_names, covs = NULL, check_order = TRUE,
                          terminal_states = c(),
                          censor_states = state_names[length(state_names)],
                          initial_states = state_names[1]) {
      checkmate::assert_class(subject_df, "data.frame")
      checkmate::assert_class(path_df, "data.frame")
      checkmate::assert_class(link_df, "data.frame")
      if (check_order) {
        path_df <- check_and_sort_paths(path_df)
      }
      if (!is.null(covs)) {
        checkmate::assert_character(covs)
      }
      if (!("draw_idx" %in% colnames(link_df))) {
        link_df$draw_idx <- 1
        link_df$draw_idx <- as.factor(link_df$draw_idx)
      }
      if (!("rep_idx" %in% colnames(link_df))) {
        link_df$rep_idx <- 1
        link_df$rep_idx <- as.factor(link_df$rep_idx)
      }
      cols1 <- c("subject_id", covs)
      cols2 <- c("path_id", "state", "time", "is_event")
      cols3 <- c("path_id", "draw_idx", "rep_idx", "subject_id")
      check_columns(subject_df, cols1)
      check_columns(path_df, cols2)
      check_columns(link_df, cols3)
      subject_df <- as_tibble(subject_df[, cols1])
      path_df <- as_tibble(path_df[, cols2])
      link_df <- as_tibble(link_df[, cols3])

      # Add the subject_index column
      ensure_numeric_factor <- function(x) {
        as.factor(as.numeric(as.factor(x)))
      }
      subject_df$subject_index <- ensure_numeric_factor(subject_df$subject_id)
      cols1 <- c("subject_index", cols1)

      # Format other
      link_df$subject_index <- ensure_numeric_factor(link_df$subject_id)
      path_df$path_id <- ensure_numeric_factor(path_df$path_id)
      link_df$path_id <- ensure_numeric_factor(link_df$path_id)
      link_df$subject_id <- NULL
      self$check_states(
        state_names = state_names, terminal_states = terminal_states,
        initial_states = initial_states, censor_states = censor_states,
        path_df = path_df
      )
      self$terminal_states <- terminal_states
      self$initial_states <- initial_states
      self$censor_states <- censor_states
      self$state_names <- state_names
      self$state_types <- seq(from = 0, to = length(state_names))
      self$covs <- covs
      self$subject_df <- subject_df
      self$path_df <- path_df
      self$link_df <- link_df
    },
    #' Get names of covariates
    covariate_names = function() {
      self$covs
    },
    #' Get path lenghts
    lengths = function() {
      self$path_df |>
        filter(is_event == TRUE) |>
        group_by(path_id) |>
        count()
    },
    #' Get longest path
    longest_path = function() {
      lens <- self$lengths()
      idx <- which(lens$n == max(lens$n))[1]
      df <- self$path_df |> filter(path_id == lens$path_id[idx])
      self$filter(unique(df$path_id))
    },
    #' Get indices of event states
    get_event_states = function() {
      match(self$get_event_state_names(), self$state_names)
    },
    get_event_state_names = function() {
      self$state_names |>
        purrr::discard(~ .x %in% self$initial_states) |>
        purrr::discard(~ .x %in% self$censor_states)
    },
    as_msdata = function(covariates = FALSE) {
      checkmate::assert_logical(covariates, len = 1)
      pathdata_to_mstate_format(self, covariates)
    },
    as_time_to_first_event = function(covs = c("subject_id"), truncate = FALSE) {
      df <- self$as_data_frame(covs, truncate = truncate)
      tt_event <- df |>
        as_time_to_first_event(states = self$get_event_states(), by = covs)
    },
    n_paths = function() {
      nrow(self$link_df)
    },
    print = function() {
      n_path <- self$n_paths()
      covs <- self$covariate_names()
      sn <- self$state_names
      message(paste0("PathData object with ", n_path, " paths"))
      message(paste0("States = {"), paste0(sn, collapse = ", "), "}")
      message(paste0("Covariates = {"), paste0(covs, collapse = ", "), "}")
    },
    check_valid_state = function(check_states, state_names, min_len) {
      checkmate::assert_character(check_states, min.len = min_len)
      assertthat::assert_that(all(check_states %in% state_names))
      assertthat::assert_that(purrr::none(check_states, duplicated))
    },
    check_states = function(state_names, terminal_states, initial_states,
                            censor_states, path_df) {
      checkmate::assert_character(state_names, min.len = 1)
      assertthat::assert_that(purrr::none(state_names, duplicated))
      observed_states <- unique(state_names[path_df$state])
      self$check_valid_state(observed_states, state_names, min_len = 1)
      self$check_valid_state(terminal_states, state_names, min_len = 0)
      self$check_valid_state(initial_states, state_names, min_len = 0)
      self$check_valid_state(censor_states, state_names, min_len = 1)
    },
    as_data_frame = function(covariates = NULL, truncate = FALSE) {
      df <- self$path_df
      if (isTRUE(truncate)) {
        term_states <- which(self$state_names %in% self$terminal_states)
        df <- df |>
          truncate_after_terminal_events(term_states)
      }
      out <- df |> inner_join(self$link_df, by = "path_id", relationship = "many-to-one")
      sub_df <- self$subject_df
      if (!is.null(covariates)) {
        sub_df <- sub_df |> select(subject_index, one_of(covariates))
      }
      out <- out |> inner_join(sub_df, by = "subject_index", relationship = "many-to-one")
      stopifnot(nrow(out) == nrow(df))
      out
    },
    as_transitions = function(only_observed = FALSE) {
      if (is.null(self$dt)) {
        self$dt <- dt <- as_transitions(self$as_data_frame(),
          state_names = self$state_names,
          state_types = self$state_types,
          terminal_states = self$terminal_states,
          censor_states = self$censor_states,
          initial_states = self$initial_states,
          covs = self$covariate_names()
        )
      } else {
        dt <- self$dt
      }
      if (only_observed) {
        # Limit possible transitions to only the observed ones
        obs <- unique(dt$df$transition)
        idx_obs <- which(dt$legend$transition %in% obs)
        n_trans <- length(idx_obs)
        dt$legend <- dt$legend[idx_obs, ]
        dt$legend$transition <- seq_len(n_trans)
        dt$df$transition <- match(
          dt$df$trans_char, dt$legend$trans_char,
          nomatch = 0
        )
      }
      dt
    },

    #' Step plot of the paths
    #'
    #' @param n_paths Number of paths to subsample for plotting.
    #' @param alpha opacity
    plot_paths = function(n_paths = NULL, alpha = 0.5) {
      df <- self$as_data_frame()
      df$is_event <- as.factor(df$is_event)
      uid <- unique(df$path_id)
      N <- length(uid)
      if (is.null(n_paths)) {
        idx_path <- seq_len(N)
      } else {
        idx_path <- sample.int(N, n_paths)
      }
      ids <- uid[idx_path]
      df <- df |>
        filter(path_id %in% ids) |>
        mutate(
          state_char = self$state_names[state],
          state = factor(state_char, levels = self$state_names, ordered = T)
        )
      ggplot(df, aes(x = time, y = state, group = path_id)) +
        geom_step(direction = "hv", alpha = alpha) +
        labs(x = "Time", y = "State", title = "State paths") +
        theme_minimal() +
        geom_point(mapping = aes(color = is_event, pch = is_event))
    },

    #' Transition proportion matrix
    #'
    #' @return a matrix
    trans_matrix = function() {
      r <- self$as_msdata()
      prop <- mstate::events(r$msdata)$Proportions
      cn <- colnames(prop)
      idx_noevent <- find_one("no event", cn)
      colnames(prop)[idx_noevent] <- "Censoring"
      prop
    },

    #' Visualize the transition proportion matrix as a graph
    #'
    #' @param digits Max number of digits to show in numbers
    #' @param ... Arguments passed to \code{\link{diagram::plotmat}}.
    #' @return a list
    plot_graph = function(digits = 3, ...) {
      f <- self$trans_matrix()
      f <- round(f, digits)
      cn <- colnames(f)
      idx_noevent <- find_one("Censoring", cn)
      idx_term <- which(colnames(f) %in% self$terminal_states)
      idx_init <- which(colnames(f) %in% self$initial_states)
      color <- rep("black", length(cn))
      col <- "gray60"
      col_term <- "firebrick"
      col_init <- "steelblue2"
      color[idx_noevent] <- col
      color[idx_term] <- col_term
      color[idx_init] <- col_init
      acol <- matrix("black", nrow(f), ncol(f))
      acol[, idx_noevent] <- col
      lcol <- acol
      lcol[, idx_term] <- col_term
      lcol[, idx_init] <- col_init
      diagram::plotmat(t(f),
        txt.col = color,
        arr.col = t(acol),
        box.lcol = t(lcol),
        lcol = t(acol),
        shadow.size = 0,
        main = "Transition proportions",
        ...
      )
    },

    # Fit CoxPH model
    coxph = function(covs = NULL, ...) {
      if (is.null(covs)) {
        covs <- self$covs
        rm <- c(
          "individual_id", "pk_post_dose", "pk_pre_dose", "country_num",
          "crcl", "weight", "dose_arm", "first_dose_amount",
          "dose_adjustment"
        )
        covs <- setdiff(covs, rm)
      }
      df <- self$as_msdata(covariates = TRUE)$msdata
      str <- paste(covs, collapse = " + ")
      form <- paste0("Surv(Tstart, Tstop, status) ~ ", str)
      survival::coxph(as.formula(form), df, ...)
    },

    # Filter based on path id, creates new object
    filter = function(path_ids_keep = NULL, subject_ids_keep = NULL,
                      rep_ids_keep = NULL, draw_ids_keep = NULL) {
      subject_df <- self$subject_df |>
        filter(is.null(subject_ids_keep) | subject_id %in% subject_ids_keep)
      link_df <- self$link_df |>
        filter(
          is.null(rep_ids_keep) | rep_idx %in% rep_ids_keep,
          is.null(!!draw_ids_keep) | draw_idx %in% draw_ids_keep,
          subject_index %in% subject_df$subject_index
        )
      path_df <- self$path_df |>
        filter(
          is.null(path_ids_keep) | path_id %in% path_ids_keep,
          path_id %in% link_df$path_id
        )
      link_df <- link_df |> filter(path_id %in% unique(path_df$path_id))
      subject_df <- subject_df |> filter(subject_index %in%
        unique(link_df$subject_index))
      link_df <- link_df |>
        left_join(subject_df |> select(subject_index, subject_id),
          by = "subject_index"
        )
      PathData$new(subject_df, path_df, link_df,
        self$state_names, self$covs,
        terminal_states = self$terminal_states,
        censor_states = self$censor_states,
        initial_states = self$initial_states
      )
    }
  )
)

shorten_name2 <- function(input_string) {
  input_string
}

as_time_to_first_event <- function(dat, states, by = c()) {
  by_syms <- rlang::syms(unique(c(by, "path_id")))
  # max time as censor time
  censor <- dat |>
    group_by(!!!by_syms) |>
    summarise(time = max(time), .groups = "keep") |>
    ungroup() |>
    expand_grid(state = states) |>
    mutate(is_event = 0)
  # min event time per state
  events <- dat |>
    filter(is_event == 1, state %in% states) |>
    group_by(state, !!!by_syms) |>
    summarize(time = min(time), .groups = "keep") |>
    ungroup() |>
    mutate(is_event = 1)
  d <- censor |>
    anti_join(events, by = c("path_id", "state")) |>
    bind_rows(events)
}

# Only creates trans_char and not actual integer index yet
as_transitions_char_single <- function(dat, state_names, terminal_states) {
  R <- nrow(dat)
  a <- dat[2:R, ]
  a$prev_state <- dat[1:(R - 1), ]$state
  s1 <- sapply(state_names[a$prev_state], shorten_name2)
  s2 <- sapply(state_names[a$state], shorten_name2)
  if (any(s1 %in% terminal_states)) {
    incorrect_states <- unique(s1[s1 %in% terminal_states])
    warning(str_c(
      "non-terminal record found for state: ", incorrect_states,
      " (which is marked as a terminal state) for subject ",
      unique(dat$subject_id), "\n"
    ))
  }
  a$trans_char <- format_transition_char(s1, s2)
  a
}

format_transition_char <- function(s1, s2) {
  paste0(s1, "->", s2)
}

# Only creates trans_char and not actual integer index yet
as_transitions_char <- function(dat, state_names, terminal_states) {
  df <- NULL
  ids <- unique(dat$subject_id)
  for (id in ids) {
    dat_id <- dat |> filter(subject_id == id)
    df <- rbind(df, as_transitions_char_single(dat_id, state_names, terminal_states))
  }

  # Edit trans_char for intervals that do not end in event
  df$trans_char[df$is_event == 0] <- NA
  df
}

# Add previous state and remove initial row of each subject
as_transitions <- function(dat, state_names, state_types, terminal_states, covs,
                           censor_states, initial_states) {
  # Transitions character representation
  dat_trans <- as_transitions_char(dat, state_names, terminal_states)

  # Create legend first
  r <- tidyr::expand_grid(prev_state = state_names, state = state_names) |>
    mutate(trans_char = format_transition_char(prev_state, state)) |>
    filter(
      !prev_state %in% terminal_states, # no transition from terminal state
      !state %in% initial_states, # no transition into initial states
      !state %in% censor_states, # no transition into or out of censor
      !prev_state %in% censor_states
    ) |>
    mutate(terminal = state %in% terminal_states) |>
    mutate(
      state = as.integer(factor(state, levels = state_names, ordered = T)),
      prev_state = as.integer(factor(prev_state, levels = state_names, ordered = T))
    ) |>
    as.data.frame()
  legend <- r[, c("trans_char", "prev_state", "state", "terminal")] |>
    filter(!is.na(trans_char))
  legend$transition <- seq_len(nrow(legend))
  trans_names <- legend$trans_char
  # add trans_type
  target_state_name <- state_names[legend$state]
  legend$trans_type <- state_types[match(target_state_name, state_names)]
  # confirm that all transitions in the data are represented in the legend
  stopifnot(all(na.omit(dat_trans$trans_char) %in% legend$trans_char))

  # Finally add actual transition index to data
  dat_trans$transition <- match(dat_trans$trans_char, trans_names)
  dat_trans$transition[which(is.na(dat_trans$transition))] <- 0

  # Order columns and drop unnecessary ones
  fields <- c(
    "time", "is_event", "state", "prev_state",
    "transition", "trans_char", "subject_id"
  )
  dat_trans <- dat_trans[, c(fields, covs)]

  # Return
  list(
    df = dat_trans,
    legend = legend
  )
}


# Function to check and sort paths based on time
check_and_sort_paths <- function(df) {
  # Initialize a flag to track if any sorting is needed
  needs_sorting <- FALSE

  # Iterate over each unique path id
  unique_paths <- unique(df$path_id)

  for (path_id in unique_paths) {
    # Get the rows for this path id
    path_rows <- df[df$path_id == path_id, ]

    # Check if the rows are sorted by time
    if (!all(order(path_rows$time) == seq_along(path_rows$time))) {
      # If not sorted, sort the rows by time
      df[df$path_id == path_id, ] <- path_rows[order(path_rows$time), ]
      # Set the flag to TRUE
      needs_sorting <- TRUE
      # Print a warning message
      warning(paste(
        "Rows for path_id", path_id,
        "were not ordered by time and have been sorted."
      ))
    }
  }

  return(df)
}

# Get subject df with numeric subject index
subject_df_with_idx <- function(pd, subs, id_map) {
  df <- pd$subject_df |> filter(subject_id %in% subs)
  df$sub_idx <- id_map$x_sub[match(df$subject_id, id_map$subject_id)]
  df
}

# Filter pathdata to subjects
filter_pathdata <- function(pd, subjects_keep) {
  df_new <- pd$df |> filter(subject_id %in% subjects_keep)
  PathData$new(df_new, pd$state_names, pd$covs,
    terminal_states = pd$terminal_states,
    censor_states = pd$censor_states,
    initial_states = pd$initial_states
  )
}

# To transition format used by mstate (msdata)
pathdata_to_mstate_format <- function(pd, covariates = FALSE) {
  dt <- pd$as_transitions()
  df <- dt$df
  siu <- unique(df$subject_id)
  df_out <- NULL
  PT <- legend_to_PT_matrix(dt$legend)
  TFI <- legend_to_TFI_matrix(dt$legend)
  for (sid in siu) {
    df_j <- df |> filter(subject_id == sid)
    df_out <- rbind(df_out, to_mstate_format_part1(df_j))
  }
  df_out <- to_mstate_format_part2(df_out, PT, TFI)
  rownames(df_out) <- NULL
  tmat <- TFI_to_mstate_transmat(TFI, pd$state_names)
  attr(df_out, "trans") <- tmat
  class(df_out) <- c("msdata", "data.frame")
  if (covariates) {
    sdf <- pd$subject_df |> mutate(id = subject_id)
    df_out <- df_out |> left_join(sdf, by = "id")
  }
  list(
    msdata = df_out,
    legend = dt$legend
  )
}

# First pass in transforming transitions format to mstate format (msdata)
# For one subject, creates the columns needed for mstate format
to_mstate_format_part1 <- function(df_sub) {
  df <- data.frame(
    id = df_sub$subject_id,
    from = df_sub$prev_state,
    to = df_sub$state,
    trans = df_sub$transition
  )
  K <- nrow(df)
  Tstart <- rep(0, K)
  if (K > 1) {
    Tstart[2:K] <- df_sub$time[1:(K - 1)]
  }
  df$Tstart <- Tstart
  df$Tstop <- df_sub$time
  df$time <- df$Tstop - df$Tstart
  df$status <- df_sub$is_event
  df
}

# Second pass in transforming transitions format to mstate format (msdata)
# Creates the additional rows corresponding to each transition
# that is at risk
to_mstate_format_part2 <- function(df, PT, TFI) {
  N <- nrow(df)
  df_out <- NULL
  for (n in 1:N) {
    row <- df[n, ]
    possible_trans <- which(as.numeric(PT[, row$from]) == 1)
    J <- length(possible_trans)
    possible_targets <- rep(0, J)
    for (j in seq_len(J)) {
      possible_targets[j] <- find_column_with_number(TFI, possible_trans[j])
    }

    if (row$status == 0) {
      rows <- NULL
      possible_targets_remaining <- possible_targets
    } else {
      rows <- row
      possible_targets_remaining <- setdiff(possible_targets, row$to)
    }

    J <- length(possible_targets_remaining)

    for (j in seq_len(J)) {
      row_rep <- row
      row_rep$to <- possible_targets_remaining[j]
      row_rep$trans <- TFI[row_rep$from, row_rep$to]
      row_rep$status <- 0
      rows <- rbind(rows, row_rep)
    }
    df_out <- rbind(df_out, rows)
  }
  df_out
}

# TFI matrix to format used by mstate
TFI_to_mstate_transmat <- function(TFI, state_names) {
  sn <- state_names[1:(length(state_names) - 1)]
  TFI[TFI == 0] <- NA
  colnames(TFI) <- sn
  rownames(TFI) <- sn
  TFI
}

fit_coxph <- function(msdat, formula_rhs = NULL) {
  if (is.null(formula_rhs)) {
    formula_rhs <- "strata(trans)"
  }
  formula <- paste0("Surv(Tstart, Tstop, status) ~ ", formula_rhs)
  coxph(as.formula(formula),
    data = msdat, method = "breslow"
  )
}
fit_mstate <- function(msdat, formula = NULL) {
  cph <- fit_coxph(msdat, formula)
  msfit(object = cph, variance = FALSE, trans = attr(msdat, "trans"))
}
fit_basehaz <- function(msdat, ...) {
  cph <- fit_coxph(msdat)
  basehaz(cph, ...)
}
plot_cumhaz_msfit <- function(msfit, legend = NULL) {
  df <- msfit$Haz
  if (!is.null(legend)) {
    leg <- legend[, c("transition", "trans_char")]
    df$transition <- df$trans
    df <- df |> left_join(leg, by = "transition")
    df$trans <- paste0(df$transition, ": ", df$trans_char)
  } else {
    df$trans <- as.factor(df$trans)
  }
  ggplot(df, aes(x = time, y = Haz, color = trans)) +
    geom_line() +
    ylab("Cumulative Hazard")
}
estimate_average_hazard <- function(msfit) {
  msfit$Haz |>
    group_by(trans) |>
    summarise(
      avg_haz = (last(Haz) - first(Haz)) / (last(time) - first(time))
    )
}

truncate_after_terminal_events <- function(df, term_states) {
  term_events <- df |>
    filter(state %in% !!term_states, is_event == 1) |>
    group_by(path_id) |>
    summarise(term_time = min(time, na.rm = T)) |>
    ungroup()
  no_terms <- df |>
    anti_join(term_events, by = "path_id")
  with_terms <- df |>
    inner_join(term_events, by = c("path_id")) |>
    filter(time <= term_time) |>
    select(-term_time)
  no_terms |>
    bind_rows(with_terms)
}

summarize_event_prob <- function(pd, target_times, by = c("subject_id")) {
  by_syms <- rlang::syms(by)
  surv_df <- as_time_to_event(pd) |>
    mutate(vv = dense_rank(str_c(!!!by_syms, sep = ":")))
  pred_pd_surv <- survival::survfit(
    survival::Surv(start_time, end_time, event) ~ vv,
    id = path_id, data = surv_df
  )

  # summarize event rates
  pp_summary <- summary(pred_pd_surv, target_times)
  pstate_cols <- pp_summary$pstate
  colnames(pstate_cols) <- pp_summary$states

  pp_cumhaz <- as_tibble(lst(
    time = pp_summary$time,
    strata = pp_summary$strata
  )) |>
    bind_cols(pstate_cols) |>
    rgeco:::.tidy_km_strata() |>
    mutate(`vv` = as.integer(vv)) |>
    left_join(surv_df |> distinct(vv, !!!by_syms), by = "vv") |>
    select(-vv)
}

# Full formula with expanded covariates
cph_full_formula <- function(msdata, ttype = FALSE) {
  a <- grep("\\.", colnames(msdata), value = TRUE)
  str <- paste(a, collapse = " + ")
  if (ttype) {
    strata <- "trans_type"
  } else {
    strata <- "trans"
  }
  paste0(str, " + strata(", strata, ")")
}



# Pathdata to single event format
to_single_event <- function(pd, event) {
  path_df <- pd$path_df
  STATE <- find_one(event, pd$state_names)
  pid <- unique(path_df$path_id)
  path_df_new <- NULL
  for (path in pid) {
    df <- path_df |>
      filter(path_id == path) |>
      arrange(time)
    df$row_num <- 1:nrow(df)
    ri <- which(df$state == STATE)
    if (length(ri) == 0) {
      df <- df[c(1, nrow(df)), ]
      df$state[2] <- 1
      df$is_event[2] <- 0
    } else {
      df <- df[c(1, ri[1]), ]
      df$state[2] <- 2
    }
    path_df_new <- rbind(path_df_new, df |> select(-row_num))
  }

  state_names <- c("Randomization", event, "Censor")
  link_df <- pd$link_df |> left_join(
    pd$subject_df |> select(subject_id, subject_index),
    by = "subject_index"
  )
  PathData$new(pd$subject_df, path_df_new, link_df,
    state_names,
    covs = pd$covs, check_order = TRUE,
    terminal_states = event
  )
}

# Look for potential covariates
potential_covariates <- function(pd, possible = NULL, ...) {
  events <- pd$get_event_state_names()
  df <- NULL
  for (e in events) {
    message("Event: ", e)
    a <- to_single_event(pd, e)
    r <- a$coxph(covs = possible, ...)
    s <- summary(r)
    pval <- summary(r)$coefficients[, 5]
    df <- rbind(df, data.frame(
      pval = as.numeric(pval), event = e,
      variable = names(pval)
    ))
  }
  df
}
