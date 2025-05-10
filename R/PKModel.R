#' One-compartment PK Model
#'
#' @export
PKModel <- R6::R6Class("PKModel",

  # PRIVATE
  private = list(
    covariates = NULL
  ),

  # PUBLIC
  public = list(

    #' @description
    #' Create model
    #'
    #' @param covariates A list with elements \code{ka},
    #' \code{CL}, and \code{V2}
    initialize = function(covariates) {
      checkmate::assert_list(covariates)
      checkmate::assert_character(covariates$ka)
      checkmate::assert_character(covariates$CL)
      checkmate::assert_character(covariates$V2)
      private$covariates <- covariates
    },

    #' @description Get the covariates affecting the k_a parameter.
    ka_covs = function() {
      private$covariates$ka
    },

    #' @description Get the covariates affecting the CL parameter.
    CL_covs = function() {
      private$covariates$CL
    },

    #' @description Get the covariates affecting the V_2 parameter.
    V2_covs = function() {
      private$covariates$V2
    },

    #' @description Format list of input PK parameters to standardized format.
    #'
    #' @param beta_pk A list of max three elements
    #' @param return A list with three elements
    format_params = function(beta_pk = NULL) {
      if (is.null(beta_pk)) {
        beta_pk <- list(ka = NULL, CL = NULL, V2 = NULL)
      }
      checkmate::assert_list(beta_pk)
      if (is.null(beta_pk$ka)) {
        x <- self$ka_covs()
        beta_pk$ka <- rep(0, length(x))
      }
      if (is.null(beta_pk$CL)) {
        x <- self$CL_covs()
        beta_pk$CL <- rep(0, length(x))
      }
      if (is.null(beta_pk$V2)) {
        x <- self$V2_covs()
        beta_pk$V2 <- rep(0, length(x))
      }
      names(beta_pk$ka) <- self$ka_covs()
      names(beta_pk$CL) <- self$CL_covs()
      names(beta_pk$V2) <- self$V2_covs()
      beta_pk
    },

    #' @description Print the object info
    #'
    #' @return nothing
    print = function() {
      x1 <- "A PKModel with:"
      x2 <- paste0(" - ka covariates: {", paste0(self$ka_covs(), collapse = ", "), "}")
      x3 <- paste0(" - CL covariates: {", paste0(self$CL_covs(), collapse = ", "), "}")
      x4 <- paste0(" - V2 covariates: {", paste0(self$V2_covs(), collapse = ", "), "}")
      msg <- paste(x1, x2, x3, x4, "\n", sep = "\n")
      cat(msg)
    },

    #' @description Simulate system in steady state
    #'
    #' @param t output time points
    #' @param theta parameter values
    #' @param dose dose
    #' @param tau dosing interval
    #' @return A vector with the same length as \code{t}, representing the
    #' concentration in the central compartment.
    simulate_ss = function(t, theta, dose, tau) {
      checkmate::assert_number(theta$ka, lower = 0)
      checkmate::assert_number(theta$CL, lower = 0)
      checkmate::assert_number(theta$V2, lower = 0)
      checkmate::assert_number(dose, lower = 0)
      checkmate::assert_number(tau, lower = 0)
      N <- length(t)
      ka <- theta$ka
      ke <- theta$CL / theta$V2
      A <- (dose / theta$V2) * (ka / (ka - ke))
      conc <- rep(0, N)
      ma <- A / (-expm1(-ka * tau))
      me <- A / (-expm1(-ke * tau))
      tt <- t %% tau
      conc <- me * exp(-ke * tt) - ma * exp(-ka * tt)
      conc
    }
  )
)

#' Compute PK params
#'
#' @param mu mean parameter values
#' @param beta values of covariate multipliers
#' @param x values of (normalized) covariates
compute_theta_pk <- function(mu, beta, x) {
  theta <- list()
  theta$ka <- exp(log(mu$ka) + sum(beta$ka * x$ka))
  theta$CL <- exp(log(mu$CL) + sum(beta$CL * x$CL))
  theta$V2 <- exp(log(mu$V2) + sum(beta$V2 * x$V2))
  theta
}
