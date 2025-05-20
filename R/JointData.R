#' Joint data class (R6 class)
#'
#' @export
#' @field paths The events data, an object of class \code{\link{PathData}}.
#' @field dosing The dosing data, an object of class \code{\link{DosingData}}.
JointData <- R6::R6Class(
  classname = "JointData",
  public = list(
    paths = NULL,
    dosing = NULL,

    #' Initialize
    #' @param paths The events data, an object of class \code{\link{PathData}}.
    #' @param dosing The dosing data, an object of class \code{\link{DosingData}}.
    initialize = function(paths, dosing) {
      checkmate::assert_class(paths, "PathData")
      checkmate::assert_class(dosing, "DosingData")
      self$paths <- paths
      self$dosing <- dosing
    },

    #' Filter subjects, creates new object
    #'
    #' @param subject_ids_keep Subjects to keep
    filter_subjects = function(subject_ids_keep) {
      checkmate::assert_character(subject_ids_keep, min.len = 1)
      new_pd <- self$paths$filter(subject_ids_keep = subject_ids_keep)
      new_dd <- self$dosing$filter(subject_ids_keep = subject_ids_keep)
      JointData$new(new_pd, new_dd)
    },

    #' @description Print info
    print = function() {
      x1 <- paste0("A JointData object:")
      cat(x1, "\n")
      self$paths$print()
      self$dosing$print()
    }
  )
)
