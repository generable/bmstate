#' Compute concordance index
#'
#' @export
#' @param df A \code{data.frame} with \code{surv} and \code{pred_prob} columns
c_index <- function(df) {
  risk_col <- "pred_prob"
  checkmate::assert_character(risk_col)
  df$p_noevent <- 1 - df[[risk_col]]
  survival::concordance(df$surv ~ df$p_noevent)
}

#' Create a data frame to evaluate state visit probability predictions
#'
#' @export
#' @param paths_obs A \code{\link{PathData}} object of observed paths.
#' @param paths_pred A \code{\link{PathData}} object of predicted paths.
#' @param state_name Name of the state of interest.
#' @param df A \code{data.frame} with \code{surv} and \code{risk} columns
create_scoring_df <- function(paths_obs, paths_pred, state_name) {
  checkmate::assert_character(state_name, len = 1)
  checkmate::assert_class(paths_obs, "PathData")
  checkmate::assert_class(paths_pred, "PathData")
  er <- p_state_visit_per_subject(paths_pred, state_name)
  a <- as_survival(paths_obs, state_name)
  a <- a |>
    dplyr::left_join(er, by = "subject_id") |>
    dplyr::rename("pred_prob" = "prob")
  a
}
