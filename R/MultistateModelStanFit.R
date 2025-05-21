#' The Fit class
#'
#' @export
#' @field stan_fit The 'Stan' fit object
#' @field stan_data Full 'Stan' data list
#' @field data A \code{\link{JointData}} object
#' @field model A \code{\link{MultistateModel}} object
MultistateModelStanFit <- R6::R6Class("MultistateModelStanFit",
  public = list(
    stan_fit = NULL,
    stan_data = NULL,
    data = NULL,
    model = NULL,

    #' @description
    #' Get model.
    #'
    #' @param name Name of model. Only has effect for
    #' \code{\link{JointModelFit}} objects.
    get_model = function(name) {
      private$model
    },

    #' @description
    #' Create model fit object
    #'
    #' @param stan_fit The 'Stan' fit object
    #' @param stan_data Full 'Stan' data list
    #' @param data A \code{\link{JointData}} object
    #' @param model A \code{\link{MultistateModel}} object
    initialize = function(model, data, stan_fit, stan_data) {
      checkmate::assert_class(data, "JointData")
      checkmate::assert_class(model, "MultistateModel")
      self$model <- model
      self$data <- data
      self$stan_fit <- stan_fit
      self$stan_data <- stan_data
    },

    #' @description Extract draws as \code{rvar}s
    #'
    #' @param name Param/quantity name
    draws = function(name = NULL) {
      d <- self$stan_fit$draws(name)
      d <- posterior::as_draws_rvars(d)
      if (is.null(name)) {
        return(d)
      }
      d[[name]]
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      x1 <- paste0("A MultistateModelStanFit.")
      msg <- paste(x1, "\n", sep = "\n")
      cat(msg)
    },

    #' @description
    #' Full names of parameters that start with \code{log_z_}.
    log_z_pars = function() {
      nams <- names(self$draws())
      match <- grepl(nams, pattern = "log_z_")
      nams[which(match)]
    },

    #' @description Generate quantities using the fit.
    #'
    #' @param stan_data Full 'Stan' input list.
    #' @param fitted_params Argument to \code{generate_quantities()}.
    #' @param ... Other arguments to \code{generate_quantities()}.
    gq = function(stan_data = NULL, fitted_params = NULL, ...) {
      if (is.null(fitted_params)) {
        fitted_params <- self$get_stan_fit()
      }
      if (inherits(fitted_params, "draws_rvars")) {
        fitted_params <- posterior::as_draws_array(fitted_params)
      }

      # Call 'Stan'
      self$get_model()$get_stanmodel()$generate_quantities(
        fitted_params = fitted_params,
        data = stan_data,
        ...
      )
    }
  )
)
