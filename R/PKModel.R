#' One-compartment PK Model
#'
#' @description Class that represents a one-compartment PK model with linear
#' oral absorption.
#' @export
PKModel <- R6::R6Class("PKModel",

  # PRIVATE
  private = list(
    covariates = NULL,
    MAX_CONC = 1e7
  ),

  # PUBLIC
  public = list(

    #' @description Get concentration upper bound
    get_max_conc = function() {
      private$MAX_CONC
    },

    #' @description Set concentration upper bound
    #' @param value Upper bound for concentration, to avoid numerical issues.
    set_max_conc = function(value) {
      checkmate::assert_number(value, lower = 0)
      message("setting max conc = ", round(value, 5))
      private$MAX_CONC <- value
      invisible(NULL)
    },

    #' @description Create model
    #'
    #' @param covariates A list with elements \code{ka},
    #' \code{CL}, and \code{V2}
    initialize = function(covariates) {
      checkmate::assert_list(covariates)
      checkmate::assert_character(covariates$ka, null.ok = TRUE)
      checkmate::assert_character(covariates$CL, null.ok = TRUE)
      checkmate::assert_character(covariates$V2, null.ok = TRUE)
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

    #' @description Get the covariates affecting the V2 parameter.
    V2_covs = function() {
      private$covariates$V2
    },

    #' @description Get names of all unique covariates
    #'
    covs = function() {
      x <- c(self$ka_covs(), self$CL_covs(), self$V2_covs())
      unique(x)
    },

    #' @description Print the object info
    #'
    #' @return nothing
    print = function() {
      x1 <- "A PKModel with:"
      x2 <- paste0(" - ka covariates: {", paste0(self$ka_covs(), collapse = ", "), "}")
      x3 <- paste0(" - CL covariates: {", paste0(self$CL_covs(), collapse = ", "), "}")
      x4 <- paste0(" - V2 covariates: {", paste0(self$V2_covs(), collapse = ", "), "}")
      MC <- self$get_max_conc()
      x5 <- paste0(" - Concentration upper bound: ", round(MC, 5))
      msg <- paste(x1, x2, x3, x4, x5, "\n", sep = "\n")
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
    },

    #' @description Compute exposure, which is the steady-state
    #' log area under concentration curve over one dosing interval
    #'
    #' @param theta parameter values
    #' @param dose dose
    #' @return A numeric value
    compute_xpsr = function(theta, dose) {
      CL <- theta[2]
      V2 <- theta[3]
      checkmate::assert_number(CL, lower = 0)
      checkmate::assert_number(V2, lower = 0)
      checkmate::assert_number(dose, lower = 0)
      log(dose) - log(CL) - log(V2)
    },

    #' @description Simulate data with many subjects
    #'
    #' @param df_subjects Data frame with one row for each subject
    #' @param beta_pk Covariate effects
    #' @param tau Dosing interval (PK)
    #' @param sigma Noise magnitude of concentration measurements (PK)
    #' @return Data frame with one row for each subject, and a
    #' \code{\link{DosingData}} object
    simulate_data = function(df_subjects, beta_pk = NULL, tau = 24,
                             sigma = 0.3) {
      checkmate::assert_class(df_subjects, "data.frame")
      checkmate::assert_true("dose" %in% colnames(df_subjects))
      checkmate::assert_number(tau, lower = 0)
      checkmate::assert_number(sigma, lower = 0)
      beta_pk <- self$format_params(beta_pk)
      N <- nrow(df_subjects)
      dd <- simulate_dosing(df_subjects, tau = tau)
      THETA <- matrix(0, N, 3)
      xu <- unique(c(self$ka_covs(), self$CL_covs(), self$V2_covs()))
      check_columns(df_subjects, xu)
      x <- data.frame(normalize_columns(as.matrix(df_subjects[, xu])))

      # Simulate observation times and parameters
      t_obs <- list()
      SUB_ID <- rep("s", N)
      for (n in seq_len(N)) {
        theta_n <- list(
          ka = exp(-3 + sum(x[n, self$ka_covs()] * beta_pk$ka) + 0.1 * stats::rnorm(1)),
          CL = exp(-2 + sum(x[n, self$CL_covs()] * beta_pk$CL) + 0.1 * stats::rnorm(1)),
          V2 = exp(-2 + sum(x[n, self$V2_covs()] * beta_pk$V2) + 0.1 * stats::rnorm(1))
        )
        t_last <- max(dd$times[[n]])
        t_pre <- t_last - (0.02 + 0.05) * runif(1) * tau
        t_post <- t_last + (0.1 + 0.2 * runif(1)) * tau
        t_obs[[n]] <- c(t_pre, t_post)
        THETA[n, ] <- unlist(theta_n)
        SUB_ID[n] <- df_subjects$subject_id[n]
      }

      # Simulate observed concentration
      CONC <- dd$simulate_pk(t_obs, THETA, self$get_max_conc())
      df_out <- NULL
      for (n in seq_len(N)) {
        xpsr <- self$compute_xpsr(THETA[n, ], dd$dose_ss[n])
        sid <- SUB_ID[n]
        conc <- (CONC |> dplyr::filter(.data$subject_id == sid))[["val"]]
        conc_noisy <- stats::rlnorm(2, meanlog = log(conc), sdlog = sigma)
        out <- c(t_obs[[n]], conc_noisy, xpsr)
        df_out <- rbind(df_out, out)
      }
      df_out <- data.frame(df_out)
      df_out <- cbind(df_out, THETA)
      colnames(df_out) <- c(
        "t_pre", "t_post", "conc_pre", "conc_post", "xpsr",
        "ka", "CL", "V2"
      )
      rownames(df_out) <- NULL
      df_out$pk_lloq <- 0
      df_out$subject_id <- df_subjects$subject_id

      # Return
      list(
        pk = df_out,
        dosing = dd,
        theta = THETA
      )
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
    }
  )
)
