#' Path data class (R6 class)
#'
#' @export
#' @description Class containing information about state paths and subjects.
#' @details It is not recommended for users to try to create data using
#' the constructor of this class. Rather use \code{\link{df_to_pathdata}}.
#' @field subject_df Data frame with one row per subject. Must have one
#' row for each subject, and
#' \code{subject_id} and all covariates as columns.
#' @field path_df Data frame of actual paths. Each row corresponds to
#' one time point. Must have \code{path_id},
#' \code{state}, \code{time}, and \code{trans_idx} as columns. These are
#' \itemize{
#'   \item \code{path_id}: path identifier
#'   \item \code{state}: integer index of the state at the time point, and
#'   until the next time point
#'   \item \code{time}: time point
#'   \item \code{trans_idx}: Integer indicating which transition this
#'   time point corresponds to. If the time point is not a state transition,
#'   this should be 0.
#' }
#' @field link_df Links the path and subject data frames. Must have
#' one row for each path, and \code{path_id}, \code{draw_idx},
#' \code{rep_idx}, and \code{subject_id} as columns.
#' @field covs Covariate column names.
#' @field transmat A \code{\link{TransitionMatrix}} describing the system to
#' which the paths belong.
PathData <- R6::R6Class(
  classname = "PathData",
  public = list(
    subject_df = NULL,
    path_df = NULL,
    link_df = NULL,
    covs = NULL,
    transmat = NULL,

    #' @description Initialize
    #' @param subject_df Data frame with one row per subject. Must have one
    #' row for each subject, and
    #' \code{subject_id} and all covariates as columns.
    #' @param path_df Data frame of actual paths. Each row corresponds to
    #' one time point. Must have \code{path_id},
    #' \code{state}, \code{time}, and \code{trans_idx} as columns. These are
    #' \itemize{
    #'   \item \code{path_id}: path identifier
    #'   \item \code{state}: integer index of the state at the time point, and
    #'   until the next time point
    #'   \item \code{time}: time point
    #'   \item \code{trans_idx}: Integer indicating which transition this
    #'   time point corresponds to. If the time point is not a state transition,
    #'   this should be 0.
    #' }
    #' @param link_df Links the path and subject data frames. Must have
    #' one row for each path, and \code{path_id}, \code{draw_idx},
    #' \code{rep_idx}, and \code{subject_id} as columns.
    #' @param covs Covariate column names.
    #' @param transmat A \code{\link{TransitionMatrix}} describing the system to
    #' which the paths belong.
    initialize = function(subject_df, path_df, link_df,
                          transmat, covs = NULL) {
      checkmate::assert_class(subject_df, "data.frame")
      checkmate::assert_class(path_df, "data.frame")
      checkmate::assert_class(link_df, "data.frame")
      checkmate::assert_class(transmat, "TransitionMatrix")
      path_df <- path_df |> dplyr::arrange(.data$path_id, .data$time)
      if (!is.null(covs)) {
        checkmate::assert_character(covs)
        covs <- setdiff(covs, "subject_id")
      }
      if (!("draw_idx" %in% colnames(link_df))) {
        link_df$draw_idx <- 1
      }
      if (!("rep_idx" %in% colnames(link_df))) {
        link_df$rep_idx <- 1
      }
      cols1 <- c("subject_id", covs)
      cols2 <- c("path_id", "state", "time", "trans_idx")
      cols3 <- c("path_id", "draw_idx", "rep_idx", "subject_id")
      check_columns(subject_df, cols1)
      check_columns(path_df, cols2)
      check_columns(link_df, cols3)
      checkmate::assert_character(subject_df$subject_id)
      checkmate::assert_character(link_df$subject_id)
      checkmate::assert_integerish(link_df$path_id)
      checkmate::assert_integerish(path_df$path_id)

      # As tibbles
      subject_df <- as_tibble(subject_df[, cols1])
      path_df <- as_tibble(path_df[, cols2])
      link_df <- as_tibble(link_df[, cols3])

      # Order
      link_df <- link_df |> dplyr::arrange(.data$path_id)
      subject_df <- subject_df |> dplyr::arrange(.data$subject_id)

      # Set fields
      self$transmat <- transmat
      self$covs <- covs
      self$subject_df <- subject_df
      self$path_df <- path_df
      self$link_df <- link_df
    },

    #' @description Get unique subject ids
    unique_subjects = function() {
      unique(self$subject_df$subject_id)
    },

    #' @description Get names of covariates
    #' @return a character vector
    covariate_names = function() {
      self$covs
    },

    #' @description Create a new PathData with a subset of covariates
    #' @param renamed_old Name of an old covariate to rename
    #' @param renamed_new Name of new covariate where old one is copied
    #' @param covs Names of all new covariates
    subset_covariates = function(covs, renamed_old = NULL, renamed_new = NULL) {
      sdf <- self$subject_df
      if (!is.null(renamed_old)) {
        checkmate::assert_character(renamed_old, len = 1)
        checkmate::assert_character(renamed_new, len = 1)
        sdf[[renamed_new]] <- sdf[[renamed_old]]
        covs <- setdiff(covs, renamed_old)
        covs <- unique(c(covs, renamed_new))
      }
      PathData$new(sdf, self$get_path_df(), self$link_df, self$transmat,
        covs = covs
      )
    },

    #' @description For each path, get the state it is in at time t
    #' @param t time
    #' @return a data frame with one row for each path
    state_at = function(t) {
      checkmate::assert_number(t, lower = 0)
      self$path_df |>
        dplyr::filter(.data$time <= t) |>
        dplyr::group_by(.data$path_id) |>
        dplyr::arrange(.data$path_id, -.data$time) |>
        dplyr::slice(1) |>
        dplyr::ungroup() |>
        dplyr::select("path_id", "state")
    },

    #' @description Get all transition times
    #' @param trans_inds Indices of transition whose occurrence times are
    #' requested. Default is \code{NULL}, in which case all transitions are
    #' considered.
    #' @return a numeric vector
    transition_times = function(trans_inds = NULL) {
      if (is.null(trans_inds)) {
        trans_inds <- self$transmat$trans_df()$trans_idx
      }
      self$as_transitions(truncate = TRUE) |>
        dplyr::filter(.data$trans_idx %in% trans_inds) |>
        dplyr::pull(.data$time)
    },

    #' @description Get number of paths
    #' @return an integer
    n_paths = function() {
      nrow(self$get_path_df() |>
        dplyr::group_by(.data$path_id) |>
        dplyr::count())
    },

    #' @description Get names of all states
    #'
    #' @return character vector
    state_names = function() {
      self$transmat$states
    },

    #' @description Get names of all states to which transitioning
    #' can be a transition event
    #'
    #' @return character vector
    get_event_state_names = function() {
      setdiff(self$state_names(), self$transmat$source_states())
    },

    #' @description Get names of terminal states
    #'
    #' @return character vector
    terminal_states = function() {
      self$transmat$absorbing_states()
    },

    #' @description Print info
    #'
    #' @return nothing
    print = function() {
      n_path <- self$n_paths()
      covs <- self$covariate_names()
      sn <- self$state_names()
      x1 <- paste0("PathData object with ", n_path, " paths")
      x2 <- paste0(" * States = {", paste0(sn, collapse = ", "), "}")
      x3 <- paste0(" * Covariates = {", paste0(covs, collapse = ", "), "}")
      msg <- paste(x1, x2, x3, "\n", sep = "\n")
      cat(msg)
    },

    #' @description Get path data frame
    #'
    #' @param truncate Remove rows after terminal events?
    get_path_df = function(truncate = FALSE) {
      checkmate::assert_logical(truncate, len = 1)
      df <- self$path_df
      if (isTRUE(truncate)) {
        term_state_idx <- self$transmat$absorbing_states(names = FALSE)
        df <- df |>
          truncate_after_terminal_events(term_state_idx)
      }
      df |> dplyr::arrange(.data$path_id, .data$time)
    },

    #' @description Convert to one long data frame
    #'
    #' @param covariates Which covariates to include?
    #' @param truncate Remove rows after terminal events?
    as_data_frame = function(covariates = NULL, truncate = FALSE) {
      fl <- self$full_link(covariates)
      df <- self$get_path_df(truncate)
      df |> dplyr::left_join(fl, by = "path_id")
    },

    #' @description Full link data frame
    #'
    #' @param covariates Which covariates to include
    #' @return A \code{data.frame} with same number of rows as \code{link_df},
    #' including also the covariate columns and \code{subject_id}
    full_link = function(covariates = NULL) {
      cols <- unique(c("subject_id", covariates))
      x <- self$subject_df[, cols]
      self$link_df[, c("subject_id", "path_id")] |>
        dplyr::left_join(x, by = "subject_id")
    },

    #' @description Data frame in transitions format
    #'
    #' @param covariates Which covariates to include?
    #' @param truncate Remove rows after terminal events first?
    #' @return A \code{data.frame}
    as_transitions = function(covariates = NULL, truncate = FALSE) {
      pdf <- self$get_path_df(truncate)
      pdf$time_prev <- c(0, pdf$time[1:(nrow(pdf) - 1)])
      tdf <- self$transmat$trans_df()
      prev_states <- c(NA, tdf$prev_state)
      states <- c(NA, tdf$state)
      out <- pdf |>
        dplyr::group_by(.data$path_id) |>
        dplyr::filter(dplyr::row_number() > 1) |>
        dplyr::ungroup() |>
        dplyr::select(
          "path_id", "time", "time_prev", "trans_idx", "state"
        )
      out$from <- prev_states[out$trans_idx + 1]
      out$to <- states[out$trans_idx + 1]
      idx_notrans <- which(out$trans_idx == 0)
      out$from[idx_notrans] <- out$state[idx_notrans]
      out$to[idx_notrans] <- out$state[idx_notrans]
      out <- out |> dplyr::select(-c("state"))
      fl <- self$full_link(covariates)
      out |> dplyr::left_join(fl, by = "path_id")
    },

    #' @description Data frame in alternative transitions format
    #'
    #' @param covariates Which covariates to include?
    #' @param truncate Remove rows after terminal events first?
    #' @return A \code{data.frame}
    as_transitions_alt = function(covariates = NULL, truncate = FALSE) {
      dt <- self$as_transitions(
        covariates = covariates,
        truncate = truncate
      )
      dt$Tstart <- dt$time_prev
      dt$Tstop <- dt$time
      dt$time <- dt$Tstop - dt$Tstart
      dt$status <- as.numeric(dt$trans_idx > 0)
      dt$trans <- dt$trans_idx
      dt |> dplyr::select(-c("time_prev", "trans_idx"))
    },

    #' @description Convert to format used by the 'mstate' package
    #'
    #' @param covariates Which covariates to include?
    #' @return An \code{msdata} object
    as_msdata = function(covariates = NULL) {
      dt <- self$as_transitions_alt(covariates, truncate = TRUE)
      df_out <- to_mstate_format(dt, self$transmat)
      attr(df_out, "trans") <- self$transmat$as_mstate_transmat()
      class(df_out) <- c("msdata", "data.frame")
      df_out
    },

    #' @description Step plot of the paths
    #'
    #' @param n_paths Number of paths to subsample for plotting.
    #' @param alpha opacity
    #' @param truncate truncate after terminal events?
    plot_paths = function(n_paths = NULL, alpha = 0.5, truncate = FALSE) {
      df <- self$as_data_frame(truncate = truncate)
      df$Transition <- as.factor(df$trans_idx > 0)
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
        labs(
          x = "Time", y = "State", title = "State paths"
        ) +
        geom_point(mapping = aes(color = .data$Transition, pch = .data$Transition))
    },

    #' @description Transition proportion matrix
    #' @param include_censor Include censoring (no event) proportion in the
    #' matrix
    #' @return a \code{table}
    prop_matrix = function(include_censor = TRUE) {
      ms <- self$as_msdata()
      prop <- mstate::events(ms)$Proportions
      if (isFALSE(include_censor)) {
        prop <- prop[, 1:(ncol(prop) - 1)]
      }
      prop
    },

    #' @description Visualize the transition proportion matrix as a graph
    #'
    #' @param digits Max number of digits to show in numbers
    #' @param ... Arguments passed to \code{qgraph}
    #' @return \code{qgraph} plot
    plot_graph = function(digits = 3, ...) {
      f <- self$prop_matrix(include_censor = FALSE)
      f <- matrix(f, ncol = ncol(f), dimnames = dimnames(f))
      transition_matrix_plot(
        f,
        self$terminal_states(),
        self$transmat$source_states(),
        edge_labs = round(f, digits = digits),
        ...
      )
    },

    #' @description Fit Cox proportional hazards model
    #'
    #' @param covariates Covariates to include.
    #' @param ... Arguments passed to \code{survival::coxph}.
    fit_coxph = function(covariates = NULL, ...) {
      msdat <- self$as_msdata(covariates = covariates)
      terms <- c("strata(trans)", covariates)
      str <- paste(terms, collapse = " + ")
      fit_coxph(msdat, formula_rhs = str, ...)
    },

    #' @description Fit frequentist 'mstate' model
    #'
    #' @param covariates Covariates to include.
    #' @param ... Arguments passed to \code{survival::coxph}.
    fit_mstate = function(covariates = NULL, ...) {
      msdat <- self$as_msdata(covariates = covariates)
      cph <- self$fit_coxph(covariates, ...)
      msdat$strata <- msdat$trans
      mstate::msfit(
        object = cph, newdata = msdat, variance = FALSE,
        trans = attr(msdat, "trans")
      )
    },

    #' @description Filter based on subject id, creates new object
    #'
    #' @param subject_ids_keep Subject ids to keep
    filter = function(subject_ids_keep) {
      checkmate::assert_character(subject_ids_keep, min.len = 1)
      subject_df <- self$subject_df |>
        dplyr::filter(subject_id %in% subject_ids_keep)
      link_df <- self$link_df |>
        dplyr::filter(subject_id %in% unique(subject_df$subject_id))
      path_df <- self$get_path_df(FALSE) |>
        dplyr::filter(path_id %in% unique(link_df$path_id))
      link_df <- link_df |> dplyr::filter(path_id %in% unique(path_df$path_id))
      subject_df <- subject_df |> dplyr::filter(subject_id %in%
        unique(link_df$subject_id))
      PathData$new(
        subject_df, path_df, link_df, self$transmat, self$covs
      )
    }
  )
)

# Fit Cox PH model using Breslow method
fit_coxph <- function(msdat, formula_rhs = NULL) {
  if (is.null(formula_rhs)) {
    formula_rhs <- "strata(trans)"
  }
  formula <- paste0("Surv(Tstart, Tstop, status) ~ ", formula_rhs)
  survival::coxph(stats::as.formula(formula),
    data = msdat, method = "breslow"
  )
}

# Plot cumulative hazard of 'msfit'
msfit_plot_cumhaz <- function(msfit, legend = NULL) {
  df <- msfit$Haz
  if (!is.null(legend)) {
    leg <- legend[, c("trans_idx", "trans_char")]
    leg$trans <- leg$trans_idx
    df$transition <- df$trans
    df <- df |> left_join(leg, by = "trans")
    df$trans <- paste0(df$transition, ": ", df$trans_char)
  } else {
    df$trans <- as.factor(df$trans)
  }
  ggplot(df, aes(x = .data$time, y = .data$Haz, color = .data$trans)) +
    geom_line() +
    ylab("Cumulative Hazard")
}

# Estimate average hazard of an 'msfit'
msfit_average_hazard <- function(msfit) {
  msfit$Haz |>
    dplyr::group_by(.data$trans) |>
    summarise(
      avg_haz = (dplyr::last(.data$Haz) - dplyr::first(.data$Haz)) /
        (dplyr::last(.data$time) - dplyr::first(.data$time))
    )
}


# Creates the additional rows corresponding to each transition
# that is at risk
to_mstate_format <- function(df, transmat) {
  N <- nrow(df)
  df_out <- NULL
  TFI <- transmat$as_transition_index_matrix()
  for (n in 1:N) {
    row <- df[n, ]
    possible_targets <- transmat$at_risk(row$from)

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


#' Look for potential covariates that affect transitions
#'
#' @export
#' @param pd A \code{\link{PathData}} object.
#' @param possible Possible covariates to look for (character vector). If
#' \code{NULL} (default) all covariates of the object are tested.
#' @param ... Arguments passed to \code{fit_coxph()}.
#' @return A \code{data.frame} with columns \code{pval}, \code{target_state},
#' \code{covariate}.
potential_covariates <- function(pd, possible = NULL, ...) {
  checkmate::assert_class(pd, "PathData")
  if (is.null(possible)) {
    possible <- pd$covariate_names()
  }
  checkmate::assert_character(possible)
  events <- pd$get_event_state_names()
  df <- NULL
  for (e in events) {
    message("Looking for covariates that affect transition to ", e)
    a <- as_single_event(pd, e)
    r <- a$fit_coxph(covariates = possible, ...)
    s <- summary(r)
    pval <- summary(r)$coefficients[, 5]
    df <- rbind(df, data.frame(
      pval = as.numeric(pval), target_state = e,
      covariate = names(pval)
    ))
  }
  df
}

#' PathData to time-to-event data format with a single event
#'
#' @export
#' @param pd A \code{\link{PathData}} object
#' @param event Name of the state corresponding to the event of interest (character)
#' @param null_state Name of the base state
#' @return A \code{\link{PathData}} object
#' @family PathData mutation functions
as_single_event <- function(pd, event, null_state = "Randomization") {
  checkmate::assert_class(pd, "PathData")
  checkmate::assert_character(event, len = 1)
  checkmate::assert_character(null_state, len = 1, min.chars = 1)
  stopifnot(event %in% pd$state_names())
  state <- which(pd$state_names() == event)
  if (length(state) != 1) {
    stop("invalid event")
  }
  path_df <- pd$get_path_df()
  pid <- unique(path_df$path_id)
  path_df_new <- NULL
  for (path in pid) {
    df <- path_df |>
      dplyr::filter(.data$path_id == path) |>
      dplyr::arrange(.data$time)
    df$row_num <- 1:nrow(df)
    ri <- which(df$state == state)
    if (length(ri) == 0) {
      df <- df[c(1, nrow(df)), ]
      df$state[2] <- 1
      df$trans_idx[2] <- 0
    } else {
      df <- df[c(1, ri[1]), ]
      df$state[2] <- 2
    }
    path_df_new <- rbind(path_df_new, df |> dplyr::select(-"row_num"))
  }
  path_df_new$trans_idx <- as.numeric(path_df_new$trans_idx > 0)

  link_df <- pd$link_df |> dplyr::left_join(
    pd$subject_df |> dplyr::select("subject_id"),
    by = "subject_id"
  )
  tm <- transmat_survival(state_names = c(null_state, event))
  PathData$new(
    pd$subject_df, path_df_new, link_df, tm, pd$covs
  )
}

#' PathData to time-to-event data format for any state other than null state
#'
#' @export
#' @inheritParams as_single_event
#' @return A \code{\link{PathData}} object
#' @family PathData mutation functions
as_any_event <- function(pd, null_state = "Randomization") {
  pd_new <- pd$clone(deep = TRUE)
  df <- pd_new$path_df
  idx <- which(pd_new$transmat$states == null_state)
  if (length(idx) != 1) {
    stop("given null_state not a state of given pd")
  }
  df$state[which(df$state == idx)] <- 1
  df$state[which(df$state != 1)] <- 2
  tm <- transmat_survival(state_names = c(null_state, "Any event"))
  pd_new$path_df <- df
  pd_new$transmat <- tm
  as_single_event(pd_new, "Any event", null_state)
}


#' PathData to event-free survival format
#'
#' @export
#' @inheritParams as_single_event
#' @return A \code{\link{PathData}} object
#' @family PathData mutation functions
as_survival <- function(pd, event) {
  checkmate::assert_class(pd, "PathData")
  N_sub <- length(pd$unique_subjects())
  a <- as_single_event(pd, event)
  ppd <- a$as_transitions()
  ppd$is_trans <- as.numeric(ppd$trans_idx > 0)
  ppd$surv <- Surv(ppd$time, ppd$is_trans)
  ppd <- ppd |> dplyr::select("path_id", "time", "surv", "is_trans", "subject_id")
  if (nrow(ppd) != N_sub) {
    stop("internal error in as_survival")
  }
  ppd
}

# Helper
count_paths_with_event <- function(c, t, S) {
  cnt <- c |>
    dplyr::filter(.data$trans_idx > 0 & .data$time <= t) |>
    dplyr::distinct(.data$path_id) |>
    dplyr::count()
  df <- data.frame(state = seq_len(S)) |> dplyr::left_join(cnt, by = "state")
  df$n[which(is.na(df$n))] <- 0
  df$n_event <- df$n
  df$n <- NULL
  df
}

#' For each subject, compute probability of visiting a given state
#' at least once before given time
#'
#' @export
#' @description Convenient wrapper for \code{\link{p_state_visit}}.
#' @inheritParams p_state_visit
#' @param state_name Name of the state (character).
#' @return A data frame
#' @family PathData summary functions
p_state_visit_per_subject <- function(pd, state_name, t = NULL) {
  checkmate::assert_character(state_name, len = 1)
  p_state_visit(pd, t, by = "subject_id") |>
    dplyr::filter(.data$state == state_name) |>
    dplyr::select(c("subject_id", "prob"))
}

#' For each non-source state, compute probability of visiting it at least once
#' before given time
#'
#' @export
#' @param pd A \code{\link{PathData}} object.
#' @param t The given time. If \code{NULL}, is set to \code{max(pd$get_path_df()$time)}.
#' @param by Factor to summarize over.
#' @return A data frame
#' @family PathData summary functions
p_state_visit <- function(pd, t = NULL, by = NULL) {
  checkmate::assert_class(pd, "PathData")
  S <- pd$transmat$num_states()
  if (is.null(t)) {
    t <- max(pd$get_path_df()$time)
  }
  checkmate::assert_numeric(t, lower = 0)

  if (!is.null(by)) {
    checkmate::assert_character(by, len = 1)
    c <- pd$as_data_frame(covariates = by) |>
      dplyr::group_by(.data$state, .data[[by]])
  } else {
    c <- pd$as_data_frame() |>
      dplyr::group_by(.data$state)
  }
  estates <- which(pd$state_names() %in% pd$get_event_state_names())
  df <- count_paths_with_event(c, t, S) |> dplyr::filter(.data$state %in% estates)
  if (!is.null(by)) {
    df_all <- c |> dplyr::ungroup()
    df_all <- df_all |> dplyr::group_by(.data[[by]])
    df_all <- df_all |>
      dplyr::distinct(.data$path_id) |>
      dplyr::count()
    df <- df |> dplyr::left_join(df_all, by = by)
    df$prob <- df$n_event / df$n
  } else {
    df$prob <- df$n_event / pd$n_paths()
  }
  df$state_idx <- df$state
  df$state <- NULL
  sdf <- pd$transmat$states_df()
  df |> dplyr::left_join(sdf, by = "state_idx")
}

# Subjects data
df_to_subjects_df <- function(dat, covs) {
  sdf <- dat[, c("subject_id", covs)]
  df_unique <- sdf |>
    dplyr::group_by(.data$subject_id) |>
    dplyr::distinct() |>
    dplyr::mutate(n_unique = dplyr::n()) |>
    dplyr::ungroup()

  if (any(df_unique$n_unique > 1)) {
    stop("Error: Some subjects have more than one unique row.")
  }

  # If no error, keep one row per subject
  df_unique |>
    dplyr::distinct(.data$subject_id, .keep_all = TRUE) |>
    dplyr::select(-"n_unique")
}

# Full df to link data frame
df_to_link_df <- function(df) {
  ldf <- df[, "subject_id"]
  ldf$path_id <- dplyr::dense_rank(ldf$subject_id)
  ldf <- ldf |>
    dplyr::group_by(.data$path_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  ldf$rep_idx <- 1L
  ldf$draw_idx <- 1L
  ldf
}

# Creates the df with path_id, state, time columns
df_to_paths_df_part1 <- function(df, link_df) {
  pdf <- df[, c("subject_id", "state", "time", "is_transition")]
  ldf <- link_df[, c("subject_id", "path_id")]
  pdf <- pdf |> dplyr::left_join(ldf, by = "subject_id")
  pdf$subject_id <- NULL
  pdf
}

# Adds the trans_idx column
df_to_paths_df_part2 <- function(pdf, tm) {
  pdf$trans_idx <- 0L
  pdf <- pdf |>
    dplyr::group_by(.data$path_id) |>
    dplyr::mutate(prev_state = dplyr::lag(.data$state, default = 0)) |>
    dplyr::ungroup()
  tim <- tm$as_transition_index_matrix()
  msg1 <- paste0(
    "Previous state should not be 0. The first row for any subject",
    " should never be marked as a transition."
  )

  for (r in seq_len(nrow(pdf))) {
    if (pdf$is_transition[r]) {
      s1 <- pdf$prev_state[r]
      if (s1 == 0) {
        stop(msg1)
      }
      s2 <- pdf$state[r]
      msg2 <- paste0(
        "row ", r, " has is_transition = TRUE but the given",
        " transition matrix has no transition from state ",
        s1, " to ", s2
      )
      t_idx <- tim[s1, s2]
      if (t_idx == 0) {
        stop(msg2)
      }
      pdf$trans_idx[r] <- t_idx
    }
  }
  pdf
}

#' Create 'PathData' from a data frame of one observed path per subject
#'
#' @export
#' @param df Data frame. Should have columns
#' \itemize{
#'   \item \code{subject_id} (character)
#'   \item \code{time} (numeric)
#'   \item \code{state} (integer, same indexing as in \code{tm})
#'   \item \code{is_transition} (logical)
#'   \item all the columns specified in \code{covs}
#' }
#' Rules:
#' \itemize{
#'    \item For each subject, there should be at least two rows, and the rows should
#'    be contiguous and the time should be non-decreasing.
#'    \item The \code{is_transition} value should indicate whether the row
#'    corresponds to a transition.
#'    \item The first row for each subject should never be a transition.
#'    \item For each row that is a transition, the transition from the state
#'    of the previous row to the current state should be a valid transition
#'    in \code{tm}.
#'    \item For each row that is not a transition, the state should not change
#'    from the previous row of that subject.
#' }
#' @param tm A \code{\link{TransitionMatrix}}
#' @param covs covariates (character vector)
#' @param validate Do stricter data validation? Recommended to use \code{TRUE}.
#' @return A \code{\link{PathData}} object
df_to_pathdata <- function(df, tm, covs = NULL, validate = TRUE) {
  check_columns(df, c("state", "time", "subject_id", "is_transition"))
  checkmate::assert_integerish(df$state)
  checkmate::assert_numeric(df$time)
  checkmate::assert_character(df$subject_id)
  checkmate::assert_logical(df$is_transition)
  df <- df |> dplyr::arrange(.data$subject_id, .data$time)
  if (validate) {
    validate_transitions(df)
  }
  checkmate::assert_class(tm, "TransitionMatrix")
  if (!is.null(covs)) {
    checkmate::assert_character(covs)
  }

  sdf <- df_to_subjects_df(df, covs)
  ldf <- df_to_link_df(df)
  pdf <- df_to_paths_df_part1(df, ldf)
  pdf <- df_to_paths_df_part2(pdf, tm) |> dplyr::select(-"prev_state")
  PathData$new(sdf, pdf, ldf, tm, covs)
}

validate_transitions <- function(df) {
  # check required columns
  if (!all(c("subject_id", "state", "is_transition") %in% names(df))) {
    stop('df must contain columns: "subject_id", "state", "is_transition"', call. = FALSE)
  }
  if (!is.logical(df$is_transition)) {
    stop('"is_transition" must be logical TRUE/FALSE.', call. = FALSE)
  }

  tmp <- df |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::group_by(.data$subject_id) |>
    dplyr::mutate(
      is_first   = dplyr::row_number() == 1,
      prev_state = dplyr::lag(.data$state),
      same_state = (.data$state == .data$prev_state) | (is.na(.data$state) & is.na(.data$prev_state))
    ) |>
    dplyr::ungroup()

  # Rule 1: first row per subject must not be a transition
  bad_first <- tmp |>
    dplyr::filter(.data$is_first & .data$is_transition)
  if (nrow(bad_first) > 0) {
    msg1 <- sprintf(
      "First row for each subject must have is_transition == FALSE (subject %s at row %d).",
      as.character(bad_first$subject_id[[1]]), bad_first$.row[[1]]
    )
    stop(msg1, call. = FALSE)
  }

  # Rule 2: non-transition rows (not first) must keep the same state
  mism <- tmp |>
    dplyr::filter(!.data$is_first & !.data$is_transition & !.data$same_state)
  if (nrow(mism) > 0) {
    msg2 <- sprintf(
      "Non-transition row changed state (subject %s at row %d).",
      as.character(mism$subject_id[[1]]), mism$.row[[1]]
    )
    stop(msg2, call. = FALSE)
  }

  invisible(TRUE)
}
