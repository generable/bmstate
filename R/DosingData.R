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
    #' @param ts A vector of output times for each subject (a list).
    #' @param theta A matrix of parameters.
    simulate_pk = function(ts, theta) {
      checkmate::assert_list(ts, len = self$num_subjects())
      pop_2cpt_partly_ss(
        ts, self$dose_ss, self$times, self$doses, theta, self$tau_ss
      )
    }
  )
)
