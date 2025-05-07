#' Create a multistate model
#'
#' @export
#' @param matrix A binary matrix where \code{matrix[i,j]} is 1 if transition
#' from state \code{i} to \code{j} is possible.
#' @param states A character vector of state names.
#' @param hazard_covs Covariates that affect the hazard. A character vector.
#' @param pk_covs Covariates that affect the PK parameters. A list with
#' elements \code{ka} \code{CL}, and \code{V2}. If \code{NULL}, a PK model
#' will not be created.
#' @param NK number of internal spline knots
#' @param compile compile 'Stan' model?
#' @return A \code{\link{MultistateModel}} object.
create_msm <- function(matrix, states, hazard_covs, pk_covs = NULL, NK = 3,
                       compile = TRUE) {
  tm <- TransitionMatrix$new(matrix, states)
  if (!is.null(pk_covs)) {
    pk <- PKModel$new(pk_covs)
  } else {
    pk <- NULL
  }

  MultistateModel$new(tm, hazard_covs, NK, pk, compile)
}

#' Main model class
#'
#' @export
#' @field transmat The transition matrix
#' @field pk_model A \code{\link{PKModel}} or NULL
MultistateModel <- R6::R6Class("MultistateModel",

  # PRIVATE
  private = list(
    stan_model = NULL,
    NK = NULL,
    hazard_covariates = NULL,
    knots = NULL,
    spline_k = 3 # degree
  ),

  # PUBLIC
  public = list(
    transmat = NULL,
    pk_model = NULL,

    #' @description
    #' Create model
    #'
    #' @param transmat A \code{\link{TransitionMatrix}}.
    #' @param P number of prediction points
    #' @param NK number of internal spline knots
    #' @param covariates The names of the hazard covariates.
    #' @param pk_model A \code{\link{PKModel}} or NULL.
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(transmat, covariates, NK, pk_model = NULL,
                          compile = TRUE) {
      checkmate::assert_integerish(NK, min = 1, len = 1)
      checkmate::assert_character(covariates)
      checkmate::assert_class(transmat, "TransitionMatrix")
      if (!is.null(pk_model)) {
        checkmate::assert_class(pk_model, "PKModel")
      }
      private$NK <- NK
      private$hazard_covariates <- covariates
      self$pk_model <- pk_model
      self$transmat <- transmat
      if (compile) {
        self$compile()
      }
    },

    #' @description Set knot locations based on event times
    #'
    #' @param t_max Max time
    #' @param t_event Occurred event times
    set_knots = function(t_max, t_event) {
      checkmate::assert_number(t_max, lower = 0)
      checkmate::assert_numeric(t_event, lower = 0, upper = t_max, min.len = 3)
      tk <- place_internal_knots(t_max, private$NK, t_event)
      private$knots <- c(0, tk, t_max)
    },

    #' @description Get knot locations
    #'
    #' @return a numeric vector
    get_knots = function() {
      if (is.null(private$knots)) {
        stop("knots have not been set, use $set_knots()")
      }
      private$knots
    },

    #' @description Get max time set for the model
    #'
    #' @return a number
    get_tmax = function() {
      max(self$get_knots())
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      covs <- self$covs()
      s <- self$transmat$states
      x1 <- paste0("A MultistateModel with:")
      x2 <- paste0(" - States: {", paste0(s, collapse = ", "), "}")
      x3 <- paste0(" - Hazard covariates: {", paste0(covs, collapse = ", "), "}")
      x4 <- paste0(" - Internal spline knots: ", private$NK)
      msg <- paste(x1, x2, x3, x4, "\n", sep = "\n")
      cat(msg)
      if (!is.null(self$pk_model)) {
        print(self$pk_model)
      }
    },

    #' @description Get the hazard covariates.
    covs = function() {
      private$hazard_covariates
    },

    #' @description Evaluate log baseline hazard
    #'
    #' @param t Output time points.
    #' @param log_w0 Intercept (log_scale).
    #' @param w Spline weights (log scale). If \code{NULL}, will be set to
    #' a vector of zeros, meaning that the log hazard is constant at
    #' \code{log_w0}
    #' @return a vector with same length as \code{t}
    log_baseline_hazard = function(t, log_w0, w = NULL) {
      tm <- self$get_tmax()
      checkmate::assert_numeric(t, lower = 0, upper = tm)
      knots <- self$get_knots()
      L <- length(knots)
      BK <- knots[c(1, L)]
      knots <- knots[2:(L - 1)]
      SBF <- bspline_basis(t, private$spline_k, knots, BK) # shape (N, L+1)
      M <- ncol(SBF)
      if (is.null(w)) {
        w <- rep(0, M)
      }
      checkmate::assert_numeric(w, len = M)
      checkmate::assert_number(log_w0)
      log_h0 <- SBF %*% w + log_w0
      as.numeric(log_h0)
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
    fit = function(pd, prior_only = FALSE, ...) {
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


# Evaluate B-spline basis at t
bspline_basis <- function(t, k, knots, BK) {
  splines2::bSpline(t,
    degree = k - 1, knots = knots, intercept = TRUE,
    Boundary.knots = BK
  )
}

# Create internal knots based on event time quantiles
place_internal_knots <- function(t_max, num_knots, t_event) {
  h <- 1 / (num_knots + 1)
  knots <- stats::quantile(t_event, probs = seq(0, 1, h))
  knots[2:(length(knots) - 1)]
}
