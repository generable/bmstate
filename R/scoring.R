#' Compute risk of event for each subject based on generated paths
#'
#' @export
#' @param p A \code{\link{PathData}} object of generated paths
#' @param event Name of the event of interest (character)
#' @return A \code{\link{PathData}} object
event_risk <- function(p, event) {
  checkmate::assert_class(p, "PathData")
  checkmate::assert_character(event, len = 1)
  stopifnot(event %in% p$get_event_state_names())
  p_event(p, by = "subject_id") |>
    dplyr::filter(.data$state == event) |>
    dplyr::select("subject_id", "prob") |>
    dplyr::rename(risk = prob)
}

#' Compute C index
#'
#' @export
#' @param df A \code{data.frame} with \code{surv} and \code{risk} columns
c_index <- function(df) {
  survival::concordance(df$surv ~ df$risk)
}
