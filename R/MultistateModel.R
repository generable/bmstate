#' Create a model
#'
#' @export
#' @param matrix a binary matrix where \code{matrix[i,j]} is 1 if transition
#' from state \code{i} to \code{j} is possible
#' @param states a character vector of state names
#' @param covariates A list with elements \code{hazard}, \code{ka},
#' \code{CL}, and \code{V2}
#' @param P number of prediction points
#' @param NK number of internal spline knots
#' @param pk should pk model be used?
#' @param compile compile 'Stan' model?
#' @return A \code{\link{MultistateModel}} object
create_msm <- function(matrix, states, covariates, P = 30, NK = 3,
                       pk = FALSE, compile = TRUE) {
  tm <- TransitionMatrix$new(matrix, states)
  MultistateModel$new(tm, covariates, P, NK, pk, compile)
}

#' Main model class
#'
#' @export
#' @field transmat The transition matrix
MultistateModel <- R6::R6Class("MultistateModel",

  # PRIVATE
  private = list(
    stan_model = NULL,
    P = NULL,
    NK = NULL,
    pk = NULL,
    covariates = NULL
  ),

  # PUBLIC
  public = list(
    transmat = NULL,

    #' @description
    #' Create model
    #'
    #' @param transmat A \code{\link{TransitionMatrix}}.
    #' @param P number of prediction points
    #' @param NK number of internal spline knots
    #' @param covariates A list with elements \code{hazard}, \code{ka},
    #' \code{CL}, and \code{V2}
    #' @param pk should pk model be used?
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(transmat, covariates, P, NK, pk, compile = TRUE) {
      checkmate::assert_integerish(P, min = 1, len = 1)
      checkmate::assert_integerish(NK, min = 1, len = 1)
      checkmate::assert_list(covariates)
      checkmate::assert_logical(pk, len = 1)
      checkmate::assert_class(transmat, "TransitionMatrix")
      private$NK <- NK
      private$P <- P
      private$covariates <- covariates
      private$pk <- pk
      self$transmat <- transmat
      if (compile) {
        self$compile()
      }
    },

    #' @description Get the covariates list.
    get_covariates = function() {
      private$covariates
    },

    #' @description Get the underlying 'Stan' model.
    get_stan_model = function() {
      private$stan_model
    },

    #' @description
    #' Create and compile the 'Stan' model.
    #'
    #' @param ... Arguments passed to \code{cmdstan_model}.
    #' @return The updated model object (invisibly).
    compile = function(...) {
      fn <- "msm.stan"
      filepath <- system.file(file.path("stan", fn), package = "bmstate")
      # silence compile warnings from cmdstan
      utils::capture.output(
        {
          mod <- cmdstanr::cmdstan_model(filepath, ...)
        },
        type = "message"
      )
      private$stan_model <- mod
      invisible(self)
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      covs <- self$get_covariates()
      message("A MultistateModel with:")
      message(" - States: {", paste0(self$transmat$states, collapse = ", "), "}")
      message(" - Hazard covariates: {", paste0(covs$hazard, collapse = ", "), "}")
      message(" - Internal spline knots: ", private$NK)
      message(" - PK: ", private$pk)
      message("   * ka covariates: {", paste0(covs$ka, collapse = ", "), "}")
      message("   * CL covariates: {", paste0(covs$CL, collapse = ", "), "}")
      message("   * V2 covariates: {", paste0(covs$V2, collapse = ", "), "}")
    },

    #' @description
    #' Create the 'Stan' data list.
    #'
    #' @param pd A \code{\link{PathData}} object.
    #' @param prior_only Sample from prior only?
    #' @return A list.
    create_standata = function(pd,
                               prior_only = FALSE) {
      checkmate::assert_class(pd, "PathData")
      checkmate::assert_logical(prior_only, len = 1)
    },

    #' @description
    #' Fit the model.
    #'
    #' @param pd A \code{\link{PathData}} object of observed paths.
    #' @param prior_only Sample from prior only.
    #' @param ... Arguments passed to \code{sample} method of the
    #' 'CmdStanR' model.
    #' @return A \code{\link{MultistateModelFit}} object.
    fit = function(pd,
                   prior_only = FALSE,
                   ...) {
      # Get Stan model object
      stan_model <- self$get_stanmodel()
      if (is.null(stan_model)) {
        stop("Stan model does not exist, you need to call compile()!")
      }

      # Create Stan input list
      d <- self$create_standata(pd, prior_only)

      # Call 'Stan'
      stan_fit <- stan_model$sample(data = d$stan_data, ...)

      # Return
      MultistateModelFit$new(self, stan_fit, pd, d$stan_data)
    }
  )
)
