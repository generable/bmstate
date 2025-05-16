#' Create the 'Stan' data list.
#'
#' @inheritParams fit_model
#' @return A list.
create_standata <- function(model, pd, prior_only = FALSE) {
  checkmate::assert_class(model, "MultistateModel")
  checkmate::assert_class(pd, "PathData")
  checkmate::assert_logical(prior_only, len = 1)
  pd
}


#' Fit the model.
#'
#' @param model A \code{\link{MultistateModel}} object.
#' @param pd A \code{\link{PathData}} object of observed paths.
#' @param prior_only Sample from prior only.
#' @param ... Arguments passed to \code{sample} method of the
#' 'CmdStanR' model.
#' @return A \code{\link{MultistateModelFit}} object.
fit_model <- function(model, pd, prior_only = FALSE, ...) {
  checkmate::assert_class(model, "MultistateModel")
  checkmate::assert_class(pd, "PathData")

  # Get Stan model object
  stan_model <- create_stan_model()

  # Create Stan input list
  d <- create_standata(model, pd, prior_only)

  # Call 'Stan'
  stan_fit <- stan_model$sample(data = d$stan_data, ...)

  # Return
  MultistateModelFit$new(self, stan_fit, pd, d$stan_data)
}
