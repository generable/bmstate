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
#' @param pd A \code{\link{PathData}} object of observed paths
#' @param risk Predicted event risk for each subject
c_index <- function(pd, risk) {
  ppd <- pd$get_path_df() |> dplyr::filter(time > 0)
  ppd$surv <- Surv(ppd$time, ppd$is_event)
  dd <- pd$link_df |> dplyr::select(path_id, subject_id)
  ppd <- ppd |>
    dplyr::left_join(dd, by = "path_id") |>
    dplyr::left_join(risk, by = "subject_id")
  ci <- survival::concordance(ppd$surv ~ ppd$risk)
  list(
    df = ppd,
    ci = ci
  )
}
