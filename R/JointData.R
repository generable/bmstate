# PathData or JointData to JointData
pd_to_jointdata <- function(data) {
  if (inherits(data, "PathData")) {
    # no dosing data, just use path data
    data <- JointData$new(data, NULL)
  }
  data
}

#' Joint data class (R6 class)
#'
#' @export
#' @description
#' Contains data about both the state paths and dosing.
#'
#' @field paths The events data, an object of class \code{\link{PathData}}.
#' @field dosing The dosing data, an object of class \code{\link{DosingData}}.
JointData <- R6::R6Class(
  classname = "JointData",
  public = list(
    paths = NULL,
    dosing = NULL,

    #' @description Initialize
    #' @param paths The events data, an object of class \code{\link{PathData}}.
    #' @param dosing The dosing data, an object of class \code{\link{DosingData}}.
    #' Can be also \code{NULL}.
    initialize = function(paths, dosing) {
      checkmate::assert_class(paths, "PathData")
      if (!is.null(dosing)) {
        checkmate::assert_class(dosing, "DosingData")
        s1 <- unique(paths$subject_df$subject_id)
        s2 <- unique(dosing$subject_ids)
        stopifnot(isTRUE(all.equal(s1, s2)))
      }
      self$paths <- paths
      self$dosing <- dosing
    },

    #' @description Filter subjects, creates new object
    #'
    #' @param subject_ids_keep Subjects to keep
    filter = function(subject_ids_keep) {
      checkmate::assert_character(subject_ids_keep, min.len = 1)
      new_pd <- self$paths$filter(subject_ids_keep = subject_ids_keep)
      if (!is.null(self$dosing)) {
        new_dd <- self$dosing$filter(subject_ids_keep = subject_ids_keep)
      } else {
        new_dd <- NULL
      }
      JointData$new(new_pd, new_dd)
    },

    #' @description Print info
    print = function() {
      x1 <- paste0("A JointData object:")
      cat(x1, "\n")
      self$paths$print()
      if (is.null(self$dosing)) {
        cat("No dosing data.\n")
      } else {
        self$dosing$print()
      }
    },

    #' @description Plot dosing
    #' @param df_fit Fit data frame.
    #' @param max_num_subjects Max number of subjects to plot.
    #' @return a \code{ggplot}
    plot_dosing = function(df_fit = NULL, max_num_subjects = 12) {
      if (is.null(self$dosing)) {
        stop("No dosing data.")
      }
      sdf <- self$paths$subject_df
      self$dosing$plot(df_fit, sdf, max_num_subjects)
    }
  )
)
