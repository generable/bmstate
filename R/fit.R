#' Fit the model using Stan.
#'
#' @param model A \code{\link{MultistateModel}} object.
#' @param data A \code{\link{JointData}} object of observed paths and dosing.
#' @param prior_only Sample from prior only?
#' @param delta_grid Time discretization delta for numerically integrating
#' hazards.
#' @param ... Arguments passed to \code{sample} method of the
#' 'CmdStanR' model.
#' @return A \code{\link{MultistateModelFit}} object.
fit_stan <- function(model, data, prior_only = FALSE,
                     delta_grid = 1, ...) {
  checkmate::assert_class(model, "MultistateModel")
  checkmate::assert_class(data, "JointData")

  # Get Stan model object
  stan_model <- create_stanmodel()

  # Create Stan input list
  d <- create_stan_data(model, data, prior_only)

  # Call 'Stan'
  stan_fit <- stan_model$sample(data = d$stan_data, ...)

  # Return
  MultistateModelFit$new(self, stan_fit, pd, d$stan_data)
}
