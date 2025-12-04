#' Dosing data class (R6 class)
#'
#' @description
#' Contains information about taken doses and dose times.
#'
#' @export
#' @field subject_ids Subject ids.
DosingData <- R6::R6Class(
  classname = "DosingData",
  public = list(
    subject_ids = NULL,

    #' @description Initialize
    #'
    #' @param subject_ids A character vector
    set_sub_ids = function(subject_ids) {
      checkmate::assert_character(subject_ids)
      self$subject_ids <- subject_ids
    },

    #' @description Print info
    print = function() {
      msg <- paste("A DosingData object with", self$num_subjects(), "subjects\n")
      cat(msg)
    },

    #' @description Get number of subjects
    num_subjects = function() {
      length(self$subject_ids)
    }
  )
)

#' Partially steady-state dosing data class (R6 class)
#'
#' @export
#' @description
#' Contains information about taken doses and dose times.
#'
#' @field doses Dose amounts (list with length equal to number of subjects).
#' Corresponds to doses after the steady state.
#' @field times Dose times (list with length equal to number of subjects).
#' Corresponds to doses after the steady state. First time here is the end
#' of steady-state assumption time range.
#' @field dose_ss Steady-state dose.
#' @field tau_ss Steady-state dosing interval.
PSSDosingData <- R6::R6Class(
  classname = "PSSDosingData",
  inherit = DosingData,
  public = list(
    doses = NULL,
    times = NULL,
    dose_ss = NULL,
    tau_ss = NULL,

    #' @description Initialize
    #'
    #' @param subject_ids A character vector
    #' @param doses Dose amounts (list with length equal to number of subjects).
    #' Corresponds to doses after the steady state.
    #' @param times Dose times (list with length equal to number of subjects).
    #' Corresponds to doses after the steady state. First time here is the end
    #' of steady-state assumption time range.
    #' @param dose_ss Steady-state dose.
    #' @param tau_ss Steady-state dosing interval.
    initialize = function(subject_ids, doses, times, dose_ss = NULL, tau_ss = 24) {
      self$set_sub_ids(subject_ids)
      N_sub <- self$num_subjects()
      checkmate::assert_number(tau_ss, lower = 0)
      checkmate::assert_list(doses, len = N_sub)
      checkmate::assert_list(times, len = N_sub)
      self$doses <- doses
      self$times <- times
      if (is.null(dose_ss)) {
        dose_ss <- rep(0, N_sub)
      }
      self$dose_ss <- dose_ss
      self$tau_ss <- tau_ss
    },

    #' @description As data frame
    #'
    as_data_frame = function() {
      N <- self$num_subjects()
      out <- NULL
      for (n in seq_len(N)) {
        sid <- self$subject_ids[n]
        dose_ss <- self$dose_ss[n]
        doses <- self$doses[[n]]
        times <- self$times[[n]]
        rows <- data.frame(
          subject_id = sid, dose_ss = dose_ss,
          dose = doses, time = times
        )
        which_missed <- which(rows$dose == 0)
        rows$dose_time <- "Taken"
        rows$dose_time[which_missed] <- "Missed"
        out <- rbind(out, rows)
      }
      out
    },

    #' @description Simulate PK dynamics
    #'
    #' @param t A vector of output times for each subject (a list).
    #' @param theta A matrix of parameters.
    #' @param skip_assert Skip most assertions and call exposed Stan directly,
    #' assuming that it exists?
    #' @param MAX_CONC concentration upper bound
    #' @return a \code{data.frame}
    simulate_pk = function(t, theta, MAX_CONC, skip_assert = FALSE) {
      checkmate::assert_list(t, len = self$num_subjects())
      checkmate::assert_logical(skip_assert, len = 1)
      checkmate::assert_number(MAX_CONC, lower = 0)
      if (skip_assert) {
        out <- pop_2cpt_partly_ss(
          t, self$dose_ss, self$times, self$doses, theta, self$tau_ss, MAX_CONC
        )
      } else {
        out <- pk_2cpt_pss(
          t, self$dose_ss, self$times, self$doses, theta, self$tau_ss, MAX_CONC
        )
      }

      time <- as.numeric(unlist(t))
      val <- unlist(out)
      val[val < 0] <- 1e-12
      sid <- rep(self$subject_ids, sapply(t, length))
      data.frame(
        time = time,
        val = val,
        subject_id = as.factor(sid)
      )
    },

    #' @description Plot dosing (and PK) data
    #'
    #' @param df_fit Fit data frame. Uses columns \code{val}, \code{lower}
    #' and \code{upper}.
    #' @param subject_df Subject data frame.
    #' @param max_num_subjects Max number of subjects to plot.
    #' @param subject_ids Which subjects to plot?
    plot = function(df_fit = NULL, subject_df = NULL, max_num_subjects = 12,
                    subject_ids = NULL) {
      dos <- self$as_data_frame()
      fitcolor <- "steelblue"
      if (is.null(max_num_subjects)) {
        max_num_subjects <- self$num_subjects()
      }
      checkmate::assert_integerish(max_num_subjects, len = 1, lower = 1)
      if (self$num_subjects() > max_num_subjects) {
        if (is.null(subject_ids)) {
          sid <- sample(unique(dos$subject_id), max_num_subjects)
        } else {
          checkmate::assert_character(subject_ids)
          sid <- subject_ids
        }

        dos <- dos |> dplyr::filter(.data$subject_id %in% sid)
        if (!is.null(df_fit)) {
          df_fit <- df_fit |> dplyr::filter(.data$subject_id %in% sid)
        }
        if (!is.null(subject_df)) {
          subject_df <- subject_df |> dplyr::filter(.data$subject_id %in% sid)
        }
      }
      plt <- ggplot(NULL, aes(
        x = .data$time, y = .data$val,
        group = .data$subject_id
      )) +
        facet_wrap(. ~ .data$subject_id, scales = "free_x") +
        geom_vline(
          data = dos, mapping = aes(xintercept = time, lty = .data$dose_time),
          col = "firebrick"
        ) +
        scale_linetype_manual(
          values =
            c(3, 1), name = "Dose"
        )
      if (!is.null(df_fit)) {
        if (!is.null(df_fit$lower)) {
          plt <- plt + geom_ribbon(
            data = df_fit, mapping = aes(ymin = lower, ymax = upper),
            fill = fitcolor, alpha = 0.7
          )
        }
        plt <- plt + geom_line(data = df_fit, color = fitcolor)
      }
      if (!is.null(subject_df)) {
        sdf <- subject_df
        df_conc <- data.frame(
          subject_id = rep(sdf$subject_id, 2),
          time = c(sdf$t_pre, sdf$t_post),
          conc = c(sdf$conc_pre, sdf$conc_post)
        )
        plt <- plt + geom_point(data = df_conc, mapping = aes(
          x = .data$time, y = .data$conc,
          group = .data$subject_id
        ))
      }

      plt + ylab("Concentration in central compartment")
    },

    #' @description Filter based on subject id, creates new object
    #'
    #' @param subject_ids_keep Subject ids to keep
    filter = function(subject_ids_keep = NULL) {
      if (is.null(subject_ids_keep)) {
        subject_ids_keep <- unique(self$subject_ids)
      }
      checkmate::assert_character(subject_ids_keep)
      idx_keep <- which(self$subject_ids %in% subject_ids_keep)
      PSSDosingData$new(
        self$subject_ids[idx_keep],
        self$doses[idx_keep],
        self$times[idx_keep],
        self$dose_ss[idx_keep],
        tau_ss = self$tau_ss
      )
    }
  )
)

#' Simulate dosing data
#'
#' @export
#' @param df_subjects Data frame with one row for each subject. Must have
#' columns \code{subject_id, num_ss_doses, num_doses, dose}.
#' @param tau Supposed dosing interval (same for each subject).
#' @param t_jitter Randomness added to dose times.
#' @param p_miss Probability of missing a dose.
#' @return A \code{\link{DosingData}} object
simulate_dosing <- function(df_subjects, tau = 24, p_miss = 0.2, t_jitter = 4) {
  checkmate::assert_data_frame(df_subjects)
  N <- nrow(df_subjects)
  needed_cols <- c("subject_id", "num_ss_doses", "num_doses", "dose")
  stopifnot(all(needed_cols %in% colnames(df_subjects)))
  doses <- list()
  times <- list()
  for (j in seq_len(N)) {
    D <- df_subjects$num_doses[j]
    Dss <- df_subjects$num_ss_doses[j]
    stopifnot(Dss < D)
    ttt <- seq(0, (D - 1) * tau, by = tau)
    ttt <- ttt + stats::rnorm(D, sd = t_jitter)
    ttt <- sort(ttt)
    ttt[1] <- 0
    ttt[which(ttt < 0)] <- 0
    if (length(unique(ttt)) != D) {
      stop("too high t_jitter?")
    }
    ddd <- rep(df_subjects$dose[j], D)
    ddd[which(stats::runif(D) < p_miss)] <- 0
    ttt <- ttt[(Dss + 1):D]
    ddd <- ddd[(Dss + 1):D]
    times[[j]] <- ttt
    doses[[j]] <- ddd
  }
  sid <- df_subjects$subject_id
  dss <- df_subjects$dose
  PSSDosingData$new(sid, doses, times, dss, tau)
}

#' Partially steady-state PK model
#'
#' @export
#' @description For each subject
#' @param t vector of output time points
#' @param dose_ss dose amount in SS (for each subject)
#' @param times time points, first of which is the end of stedy-state assumption
#' @param doses doses taken after \code{t_last_ss}
#' @param theta PK params for each subject
#' @param tau Dosing interval (same for all subjects).
#' @param MAX_CONC Concentration upper bound.
#' @return For each subject, the concentration in the central compartment at
#' times \code{t}
pk_2cpt_pss <- function(t, dose_ss, times, doses, theta, tau, MAX_CONC) {
  ensure_exposed_stan_functions()
  checkmate::assert_number(tau, lower = 0)
  checkmate::assert_number(MAX_CONC, lower = 0)
  checkmate::assert_numeric(dose_ss, lower = 0)
  N_sub <- length(dose_ss)
  checkmate::assert_list(times, len = N_sub)
  checkmate::assert_list(doses, len = N_sub)
  checkmate::assert_list(t, len = N_sub)
  checkmate::assert_matrix(theta, nrows = N_sub, ncols = 3)
  pop_2cpt_partly_ss(t, dose_ss, times, doses, theta, tau, MAX_CONC)
}
