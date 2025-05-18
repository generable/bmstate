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
      data.frame(time = time, val = val, subject_id = as.factor(sid))
    },

    #' Plot PK data simulated using the dosing schedule
    #'
    #' @param df A data frame simulated using \code{simulate_pk}.
    plot_pk = function(df) {
      plt <- ggplot(df, aes(x = .data$time, y = .data$val, group = .data$subject_id)) +
        facet_wrap(. ~ .data$subject_id) +
        geom_line()
      plt + ylab("Concentration in central compartment")
    }
  )
)
