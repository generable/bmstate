#' Dosing data class (R6 class)
#'
#' @export
#' @field doses Dose amounts (list with length equal to number of subjects).
#' @field times Dose times (list with length equal to number of subjects).
#' @field dose_ss Steady-state dose.
#' @field tau_ss Steady-state dosing interval.
#' @field subject_ids Subject ids.
DosingData <- R6::R6Class(
  classname = "DosingData",
  public = list(
    doses = NULL,
    times = NULL,
    dose_ss = NULL,
    tau_ss = NULL,
    subject_ids = NULL,

    #' Initialize
    #'
    #' @param subject_ids A character vector
    #' @param doses Dose amounts (list with length equal to number of subjects).
    #' @param times Dose times (list with length equal to number of subjects).
    #' @param dose_ss Steady-state dose.
    #' @param tau_ss Steady-state dosing interval.
    initialize = function(subject_ids, doses, times, dose_ss = NULL, tau_ss = 24) {
      checkmate::assert_character(subject_ids)
      N_sub <- length(subject_ids)
      checkmate::assert_number(tau_ss, lower = 0)
      checkmate::assert_list(doses, len = N_sub)
      N_sub <- length(doses)
      checkmate::assert_list(times, len = N_sub)
      self$subject_ids <- subject_ids
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
        out <- rbind(out, rows)
      }
      out
    },

    #' @description Print info
    print = function() {
      msg <- paste("A DosingData object with", self$num_subjects(), "subjects\n")
      cat(msg)
    },

    #' @description Get number of subjects
    num_subjects = function() {
      length(self$subject_ids)
    },

    #' Simulate PK dynamics
    #'
    #' @param t A vector of output times for each subject (a list).
    #' @param theta A matrix of parameters.
    #' @return a \code{data.frame}
    simulate_pk = function(t, theta) {
      checkmate::assert_list(t, len = self$num_subjects())
      out <- pop_2cpt_partly_ss(
        t, self$dose_ss, self$times, self$doses, theta, self$tau_ss
      )
      time <- as.numeric(unlist(t))
      val <- unlist(out)
      sid <- rep(self$subject_ids, sapply(t, length))
      data.frame(
        time = time,
        val = val,
        subject_id = as.factor(sid)
      )
    },

    #' Plot dosing (and PK) data
    #'
    #' @param df_fit Fit data frame.
    #' @param subject_df Subject data frame.
    #' @param max_num_subjects Max number of subjects to plot.
    plot = function(df_fit = NULL, subject_df = NULL, max_num_subjects = 20) {
      dos <- self$as_data_frame()
      fitcolor <- "steelblue"
      if (is.null(max_num_subjects)) {
        max_num_subjects <- self$num_subjects()
      }
      checkmate::assert_integerish(max_num_subjects, len = 1, lower = 1)
      if (self$num_subjects() > max_num_subjects) {
        sid <- sample(unique(dos$subject_id), max_num_subjects)
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
        facet_wrap(. ~ .data$subject_id) +
        geom_vline(
          data = dos, mapping = aes(xintercept = time),
          col = "firebrick", lty = 2
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
      DosingData$new(
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
#' @param df_subjects Data frame with one row for each subject
#' @param tau Dosing interval.
#' @return A \code{\link{DosingData}} object
simulate_dosing <- function(df_subjects, tau = 24) {
  N <- nrow(df_subjects)
  dose_ss <- c(15, 30, 60)[sample.int(3, N, replace = TRUE)]
  t1 <- 100 + 100 * stats::runif(N)
  t2 <- t1 + (1 - 0.5 * stats::runif(N)) * tau
  d1 <- c(30, 0)[sample.int(2, N, replace = TRUE)]
  d2 <- c(60, 0)[sample.int(2, N, replace = TRUE)]
  times <- as.list(data.frame(t(matrix(c(t1, t2), ncol = 2))))
  doses <- as.list(data.frame(t(matrix(c(d1, d2), ncol = 2))))
  sid <- df_subjects$subject_id
  DosingData$new(sid, doses, times, dose_ss, tau)
}
