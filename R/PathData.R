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

    #' @description Get path lengths (among paths that include events)
    #' @param truncate Remove rows after terminal events first?
    #' @return a data frame of path ids and counts
    lengths = function(truncate = FALSE) {
      self$get_path_df(truncate) |>
        dplyr::filter(is_event == TRUE) |>
        dplyr::group_by(path_id) |>
        dplyr::count()
    },

    #' @description Get number of paths
    #' @return an integer
    n_paths = function() {
      nrow(self$get_path_df() |>
        dplyr::group_by(.data$path_id) |>
        dplyr::count())
    },

    #' @description Get longest path
    #'
    #' @return a \code{\link{PathData}} object with just one path
    longest_path = function() {
      lens <- self$lengths(FALSE)
      idx <- which(lens$n == max(lens$n))[1]
      df <- self$get_path_df(FALSE) |> dplyr::filter(path_id == lens$path_id[idx])
      self$filter(unique(df$path_id))
    },

    #' @description Get name of null state
    null_state = function() {
      self$transmat$source_states()
    },

    #' @description Get names of all states
    #'
    #' @return character vector
    state_names = function() {
      self$transmat$states
    },

    #' @description Get names of terminal states
    #'
    #' @return character vector
    terminal_states = function() {
      self$transmat$absorbing_states()
    },

    #' @description Get indices of event states
    #'
    #' @return integer vector
    get_event_states = function() {
      match(self$get_event_state_names(), self$state_names())
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
      setdiff(self$state_names(), self$null_state())
    },

    #' @description Print info
    #'
    #' @return nothing
    print = function() {
      n_path <- self$n_paths()
      covs <- self$covariate_names()
      sn <- self$get_event_state_names()
      x1 <- paste0("PathData object with ", n_path, " paths")
      x2 <- paste0(" * Null state = {", self$null_state(), "}")
      x3 <- paste0(" * Event states = {", paste0(sn, collapse = ", "), "}")
      x4 <- paste0(" * Covariates = {", paste0(covs, collapse = ", "), "}")
      msg <- paste(x1, x2, x3, x4, "\n", sep = "\n")
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
    #' including also the covariate columns and \code{subject_index}
    full_link = function(covariates = NULL) {
      x <- self$subject_df[, c("subject_index", covariates)]
      self$link_df[, c("subject_index", "path_id")] |>
        dplyr::left_join(x, by = "subject_index")
    },

    #' @description Data frame in transitions format
    #'
    #' @param covariates Which covariates to include?
    #' @param truncate Remove rows after terminal events first?
    #' @return A \code{data.frame}
    as_transitions = function(covariates = NULL, truncate = FALSE) {
      pdf <- self$get_path_df(truncate)
      pdf$time_prev <- c(0, pdf$time[1:(nrow(pdf) - 1)])
      pdf$keep <- as.numeric((pdf$trans_idx > 0) | (pdf$is_censor == 1))
      tdf <- self$transmat$trans_df()
      prev_states <- c(NA, tdf$prev_state)
      states <- c(NA, tdf$state)
      out <- pdf |>
        dplyr::filter(.data$keep == 1) |>
        dplyr::select(
          "path_id", "time", "time_prev", "trans_idx", "is_censor",
          "state", "is_event"
        )
      out$from <- prev_states[out$trans_idx + 1]
      out$to <- states[out$trans_idx + 1]
      idx_censor <- which(out$is_censor == 1)
      out$from[idx_censor] <- out$state[idx_censor]
      out$to[idx_censor] <- out$state[idx_censor]
      out <- out |> dplyr::select(-c("state", "is_censor"))
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
      dt$status <- dt$is_event
      dt$trans <- dt$trans_idx
      dt |> dplyr::select(-c("is_event", "time_prev", "trans_idx"))
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
        self$null_state(),
        edge_labs = TRUE,
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
      mstate::msfit(object = cph, variance = FALSE, trans = attr(msdat, "trans"))
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
      path_df <- self$get_path_df(FALSE) |>
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
      PathData$new(
        subject_df, path_df, link_df,
        self$transmat,
        self$covs
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

#' Fit Cox PH model using Breslow method
#'
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

#' Plot cumulative hazard of 'msfit'
#'
#' @export
#' @param msfit An \code{msfit} object
#' @param legend transition name legend
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

#' Estimate average hazard of an 'msfit'
#'
#' @export
#' @param msfit An \code{msfit} object
msfit_average_hazard <- function(msfit) {
  msfit$Haz |>
    dplyr::group_by(trans) |>
    summarise(
      avg_haz = (dplyr::last(Haz) - dplyr::first(Haz)) /
        (dplyr::last(time) - dplyr::first(time))
    )
}


# Creates the additional rows corresponding to each transition
# that is at risk
to_mstate_format <- function(df, transmat) {
  df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      possible_targets = list(transmat$at_risk(from)),
      remaining_targets = list(setdiff(possible_targets, to)),
      all_targets = list(c(to, remaining_targets)),
      status_all = list(c(status, rep(0, length(remaining_targets))))
    ) |>
    dplyr::select(-"to", -"status", -"possible_targets", -"remaining_targets") |>
    tidyr::unnest(c(all_targets, status_all)) |>
    dplyr::rename(to = all_targets, status = status_all) |>
    dplyr::ungroup()
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
