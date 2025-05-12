# util
check_columns <- function(df, needed_columns) {
  if (!(all(needed_columns %in% colnames(df)))) {
    message("found = {", paste(colnames(df), collapse = ", "), "}")
    message("required = {", paste(needed_columns, collapse = ", "), "}")
    stop("some needed columns are missing from df")
  }
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
#' @field transmat A \code{\link{TransitionMatrix}} describing the system in to
#' which the paths belong.
PathData <- R6::R6Class(
  classname = "PathData",
  private = list(
    dt = NULL
  ),
  public = list(
    subject_df = NULL,
    path_df = NULL,
    link_df = NULL,
    covs = NULL,
    transmat = NULL,

    #' Initialize
    #' @param subject_df Data frame with one row per subject. Must have
    #' \code{subject_id} and all covariates as columns.
    #' @param path_df Data frame of actual paths. Must have \code{path_id},
    #' \code{state}, \code{time}, \code{is_event}, \code{is_censor}
    #' and \code{trans_idx} as columns.
    #' @param link_df Links the path and subject data frames. Must have
    #' \code{path_id}, \code{draw_idx}, \code{rep_idx}, and \code{subject_id} as
    #' columns.
    #' @param covs Covariate column names.
    #' @param check_order Check order of paths?
    #' @param transmat A \code{\link{TransitionMatrix}} describing the system to
    #' which the paths belong.
    initialize = function(subject_df, path_df, link_df,
                          transmat, covs = NULL, check_order = TRUE) {
      checkmate::assert_class(subject_df, "data.frame")
      checkmate::assert_class(path_df, "data.frame")
      checkmate::assert_class(link_df, "data.frame")
      checkmate::assert_class(transmat, "TransitionMatrix")
      if (check_order) {
        path_df <- check_and_sort_paths(path_df)
      }
      if (!is.null(covs)) {
        checkmate::assert_character(covs)
        covs <- setdiff(covs, "subject_id")
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
      cols2 <- c("path_id", "state", "time", "is_event", "is_censor", "trans_idx")
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
      self$transmat <- transmat
      self$covs <- covs
      self$subject_df <- subject_df
      self$path_df <- path_df
      self$link_df <- link_df
    },

    #' @description Get names of covariates
    #' @return a character vector
    covariate_names = function() {
      self$covs
    },
    #' @description Get path lenghts (among paths that include events)
    #' @return a data frame of path ids and counts
    lengths = function() {
      self$path_df |>
        dplyr::filter(is_event == TRUE) |>
        dplyr::group_by(path_id) |>
        dplyr::count()
    },
    #' @description Get number of paths
    #' @return an integer
    n_paths = function() {
      nrow(self$path_df |>
        dplyr::group_by(.data$path_id) |>
        dplyr::count())
    },
    #' @description Get longest path
    #' @return a \code{\link{PathData}} object with just one path
    longest_path = function() {
      lens <- self$lengths()
      idx <- which(lens$n == max(lens$n))[1]
      df <- self$path_df |> dplyr::filter(path_id == lens$path_id[idx])
      self$filter(unique(df$path_id))
    },
    #' @description Get name of censoring state
    censor_state = function() {
      self$transmat$censor_state
    },
    #' @description Get name of null state
    null_state = function() {
      self$transmat$source_states()
    },
    #' @description Get names of all states, including censor
    #' @return character vector
    state_names = function() {
      c(self$transmat$states, self$censor_state())
    },
    #' @description Get names of terminal states
    #' @return character vector
    terminal_states = function() {
      self$transmat$terminal_states()
    },
    #' @description Get indices of event states
    #' @return integer vector
    get_event_states = function() {
      match(self$get_event_state_names(), self$state_names())
    },

    #' @description Convert to format used by the 'mstate' package
    #'
    #' @param covariates Include covariates?
    #' @return A list with elements \code{msdata} and
    #' \code{legend}
    as_msdata = function(covariates = FALSE) {
      checkmate::assert_logical(covariates, len = 1)
      pathdata_to_mstate_format(self, covariates)
    },
    #' @description Format as time to first event
    #'
    #' @param covs covariates to include in output data frame
    #' @param truncate truncate after terminal events?
    as_time_to_first_event = function(covs = c("subject_id"), truncate = FALSE) {
      df <- self$as_data_frame(covs, truncate = truncate)
      df |>
        as_time_to_first_event(states = self$get_event_states(), by = covs)
    },

    #' @description Get names of event states
    #'
    #' @return a character vector
    get_event_state_names = function() {
      setdiff(self$state_names(), c(self$censor_state(), self$null_state()))
    },

    #' @description Print info
    #'
    #' @return nothing
    print = function() {
      n_path <- self$n_paths()
      covs <- self$covariate_names()
      sn <- self$get_event_state_names()

      message(paste0("PathData object with ", n_path, " paths"))
      message(paste0(" * Null state = {", self$null_state(), "}"))
      message(paste0(" * Censoring state = {", self$censor_state(), "}"))
      message(paste0(" * Event states = {"), paste0(sn, collapse = ", "), "}")
      message(paste0(" * Covariates = {"), paste0(covs, collapse = ", "), "}")
    },

    #' @description Convert to one long data frame
    #'
    #' @param covariates Which covariates to include?
    #' @param truncate Truncate after terminal events?
    as_data_frame = function(covariates = NULL, truncate = FALSE) {
      df <- self$path_df
      if (isTRUE(truncate)) {
        term_state_idx <- self$transmat$absorbing_states(names = FALSE)
        df <- df |>
          truncate_after_terminal_events(term_state_idx)
      }
      out <- df |>
        inner_join(self$link_df, by = "path_id", relationship = "many-to-one")
      sub_df <- self$subject_df
      if (!is.null(covariates)) {
        sub_df <- sub_df |> dplyr::select(subject_index, one_of(covariates))
      }
      out <- out |>
        inner_join(sub_df, by = "subject_index", relationship = "many-to-one")
      stopifnot(nrow(out) == nrow(df))
      out
    },

    #' @description Data frame in transitions format
    #'
    #' @param only_observed Limit transitions to only the observed ones?
    #' @param force_rerun Force rerun of the conversion? If \code{FALSE},
    #' a cached version is used if it exists.
    #' @return A list with elements \code{df} (transition data frame) and
    #' \code{legend}
    as_transitions = function(only_observed = FALSE, force_rerun = FALSE) {
      if (is.null(private$dt) || force_rerun) {
        private$dt <- dt <- as_transitions(self$as_data_frame(),
          state_names = self$state_names,
          state_types = seq(from = 0, to = length(self$state_names)),
          terminal_states = self$terminal_states,
          censor_state = self$censor_state,
          null_state = self$null_state,
          covs = self$covariate_names()
        )
      } else {
        dt <- private$dt
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

    #' @description Step plot of the paths
    #'
    #' @param n_paths Number of paths to subsample for plotting.
    #' @param alpha opacity
    #' @param truncate truncate after terminal events?
    plot_paths = function(n_paths = NULL, alpha = 0.5, truncate = FALSE) {
      df <- self$as_data_frame(truncate = truncate)
      df$is_event <- as.factor(df$is_event)
      sn <- self$state_names()
      uid <- unique(df$path_id)
      N <- length(uid)
      if (is.null(n_paths)) {
        idx_path <- seq_len(N)
      } else {
        idx_path <- sample.int(N, n_paths)
      }
      ids <- uid[idx_path]
      df <- df |>
        dplyr::filter(.data$path_id %in% ids) |>
        mutate(
          state_char = sn[state],
          state = factor(state_char, levels = sn, ordered = T)
        )
      ggplot(df, aes(x = .data$time, y = .data$state, group = .data$path_id)) +
        geom_step(direction = "hv", alpha = alpha) +
        labs(x = "Time", y = "State", title = "State paths") +
        geom_point(mapping = aes(color = .data$is_event, pch = .data$is_event))
    },

    #' @description Transition proportion matrix
    #'
    #' @return a matrix
    trans_matrix = function() {
      r <- self$as_msdata()
      prop <- mstate::events(r$msdata)$Proportions
      cn <- colnames(prop)
      idx_noevent <- find_one("no event", cn)
      colnames(prop)[idx_noevent] <- self$censor_state
      prop <- rbind(prop, rep(0, ncol(prop)))
      rownames(prop)[nrow(prop)] <- self$censor_state
      prop
    },

    #' @description Visualize the transition proportion matrix as a graph
    #'
    #' @param digits Max number of digits to show in numbers
    #' @param ... Arguments passed to \code{qgraph}
    #' @param include_censor Include censoring state in the graph?
    #' @return \code{qgraph} plot
    plot_graph = function(digits = 3, include_censor = FALSE, ...) {
      f <- self$trans_matrix()
      transition_matrix_plot(
        f, self$terminal_states,
        self$censor_state, self$null_state,
        include_censor,
        edge_labs = TRUE, ...
      )
    },

    #' @description Fit Cox proportional hazards model
    #'
    #' @param covs Covariate names.
    #' @param ... Arguments passed to \code{survival::coxph}.
    coxph = function(covs = NULL, ...) {
      if (is.null(covs)) {
        covs <- self$covs
        rm <- c(
          "individual_id", "pk_post_dose", "pk_pre_dose", "country_num",
          "crcl", "weight", "dose_arm", "first_dose_amount",
          "dose_adjustransmatent"
        )
        covs <- setdiff(covs, rm)
      }
      df <- self$as_msdata(covariates = TRUE)$msdata
      str <- paste(covs, collapse = " + ")
      form <- paste0("Surv(Tstart, Tstop, status) ~ ", str)
      survival::coxph(stats::as.formula(form), df, ...)
    },

    #' @description Filter based on path id, creates new object
    #'
    #' @param path_ids_keep Path ids to keep
    #' @param subject_ids_keep Subject ids to keep
    #' @param rep_ids_keep Repetition ids to keep
    #' @param draw_ids_keep Draw ids to keep
    filter = function(path_ids_keep = NULL, subject_ids_keep = NULL,
                      rep_ids_keep = NULL, draw_ids_keep = NULL) {
      subject_df <- self$subject_df |>
        dplyr::filter(is.null(subject_ids_keep) | subject_id %in% subject_ids_keep)
      link_df <- self$link_df |>
        dplyr::filter(
          is.null(rep_ids_keep) | rep_idx %in% rep_ids_keep,
          is.null(!!draw_ids_keep) | draw_idx %in% draw_ids_keep,
          subject_index %in% subject_df$subject_index
        )
      path_df <- self$path_df |>
        dplyr::filter(
          is.null(path_ids_keep) | path_id %in% path_ids_keep,
          path_id %in% link_df$path_id
        )
      link_df <- link_df |> dplyr::filter(path_id %in% unique(path_df$path_id))
      subject_df <- subject_df |> dplyr::filter(subject_index %in%
        unique(link_df$subject_index))
      link_df <- link_df |>
        left_join(subject_df |> dplyr::select(subject_index, subject_id),
          by = "subject_index"
        )
      PathData$new(subject_df, path_df, link_df,
        self$state_names, self$covs,
        terminal_states = self$terminal_states,
        censor_state = self$censor_state,
        null_state = self$null_state
      )
    }
  )
)

as_time_to_first_event <- function(dat, states, by = c()) {
  by_syms <- rlang::syms(unique(c(by, "path_id")))
  # max time as censor time
  censor <- dat |>
    dplyr::group_by(!!!by_syms) |>
    summarise(time = max(time), .groups = "keep") |>
    dplyr::ungroup() |>
    expand_grid(state = states) |>
    mutate(is_event = 0)
  # min event time per state
  events <- dat |>
    dplyr::filter(is_event == 1, state %in% states) |>
    dplyr::group_by(state, !!!by_syms) |>
    summarize(time = min(time), .groups = "keep") |>
    dplyr::ungroup() |>
    mutate(is_event = 1)
  censor |>
    anti_join(events, by = c("path_id", "state")) |>
    dplyr::bind_rows(events)
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

#' Get subject df with numeric subject index
#'
#' @export
#' @param pd A \code{\link{PathData}} object
#' @param subs Subject indices (from \code{\link{do_split}})
#' @param id_map Maps numeric id to original id
subject_df_with_idx <- function(pd, subs, id_map) {
  checkmate::assert_class(pd, "PathData")
  df <- pd$subject_df |> dplyr::filter(subject_id %in% subs)
  df$sub_idx <- id_map$x_sub[match(df$subject_id, id_map$subject_id)]
  df
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
    df_j <- df |> dplyr::filter(subject_id == sid)
    df_out <- rbind(df_out, to_mstate_format_part1(df_j))
  }
  df_out <- to_mstate_format_part2(df_out, PT, TFI)
  rownames(df_out) <- NULL
  transmatat <- TFI_to_mstate_transmat(TFI, pd$state_names)
  attr(df_out, "trans") <- transmatat
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

#' Fit Cox PH model using Breslow method
#'
#' @export
#' @param msdat \code{msdata} object
#' @param formula_rhs Formula right hand side that is appended to
#' \code{Surv(Tstart, Tstop, status) ~ }. If \code{NULL} (default), then
#' \code{strata(trans)} is used
#' @return value returned by \code{survival::coxph}
fit_coxph <- function(msdat, formula_rhs = NULL) {
  if (is.null(formula_rhs)) {
    formula_rhs <- "strata(trans)"
  }
  formula <- paste0("Surv(Tstart, Tstop, status) ~ ", formula_rhs)
  survival::coxph(stats::as.formula(formula),
    data = msdat, method = "breslow"
  )
}

#' Fit 'mstate' model
#'
#' @export
#' @param msdat \code{msdata} object
#' @param formula formula
fit_mstate <- function(msdat, formula = NULL) {
  cph <- fit_coxph(msdat, formula)
  mstate::msfit(object = cph, variance = FALSE, trans = attr(msdat, "trans"))
}
fit_basehaz <- function(msdat, ...) {
  cph <- fit_coxph(msdat)
  survival::basehaz(cph, ...)
}

#' Plot cumulative hazard of 'msfit'
#'
#' @export
#' @param msfit An \code{msfit} object
#' @param legend transition name legend
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

#' Estimate average hazard of an 'msfit'
#'
#' @export
#' @param msfit An \code{msfit} object
estimate_average_hazard <- function(msfit) {
  msfit$Haz |>
    dplyr::group_by(trans) |>
    summarise(
      avg_haz = (dplyr::last(Haz) - dplyr::first(Haz)) /
        (dplyr::last(time) - dplyr::first(time))
    )
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
    dplyr::bind_cols(pstate_cols) |>
    rgeco:::.tidy_km_strata() |>
    mutate(`vv` = as.integer(vv)) |>
    left_join(surv_df |> dplyr::distinct(vv, !!!by_syms), by = "vv") |>
    dplyr::select(-vv)
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
      dplyr::filter(path_id == path) |>
      dplyr::arrange(time)
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
    path_df_new <- rbind(path_df_new, df |> dplyr::select(-row_num))
  }

  state_names <- c("Randomization", event, "Censor")
  link_df <- pd$link_df |> left_join(
    pd$subject_df |> dplyr::select(subject_id, subject_index),
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
