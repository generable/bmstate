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
    },

    #' @description Compute steady-state area under curve over one dosing interval
    #'
    #' @param theta parameter values
    #' @param dose dose
    #' @return A numeric value
    compute_ss_auc = function(theta, dose) {
      CL <- theta[2]
      checkmate::assert_number(CL, lower = 0)
      checkmate::assert_number(dose, lower = 0)
      dose / CL
    },

    #' @description Simulate data with many subjects
    #'
    #' @param df_subjects Data frame with one row for each subject
    #' @param beta_pk Covariate effects
    #' @param tau Dosing interval
    #' @param sigma Noise magnitude
    #' @return Data frame with one row for each subject, and a
    #' \code{\link{DosingData}} object
    simulate_data = function(df_subjects, beta_pk = NULL, tau = 24,
                             sigma = 0.2) {
      checkmate::assert_class(df_subjects, "data.frame")
      checkmate::assert_true("dose" %in% colnames(df_subjects))
      checkmate::assert_number(tau, lower = 0)
      checkmate::assert_number(sigma, lower = 0)
      beta_pk <- self$format_params(beta_pk)
      N <- nrow(df_subjects)
      dd <- simulate_dosing(df_subjects, tau = tau)
      THETA <- matrix(0, N, 3)

      # Simulate observation times and parameters
      t_obs <- list()
      SUB_ID <- rep("s", N)
      for (n in seq_len(N)) {
        row <- df_subjects[n, ]
        theta_n <- list(
          ka = exp(-2 + sum(row[, self$ka_covs()] * beta_pk$ka)),
          CL = exp(-2 + sum(row[, self$CL_covs()] * beta_pk$CL)),
          V2 = exp(-2 + sum(row[, self$V2_covs()] * beta_pk$V2))
        )
        t_last <- max(dd$times[[n]])
        t_pre <- t_last - (0.02 + 0.05) * runif(1) * tau
        t_post <- t_last + (0.1 + 0.2 * runif(1)) * tau
        t_obs[[n]] <- c(t_pre, t_post)
        THETA[n, ] <- unlist(theta_n)
        SUB_ID[n] <- row$subject_id
      }

      # Simulate observed concentration
      CONC <- dd$simulate_pk(t_obs, THETA)
      df_out <- NULL
      for (n in seq_len(N)) {
        ss_auc <- self$compute_ss_auc(THETA[n, ], dd$dose_ss[n])
        sid <- SUB_ID[n]
        conc <- (CONC |> dplyr::filter(.data$subject_id == sid))[["val"]]
        conc_noisy <- stats::rlnorm(2, meanlog = log(conc), sdlog = sigma)
        out <- c(t_obs[[n]], conc_noisy, ss_auc)
        df_out <- rbind(df_out, out)
      }
      df_out <- data.frame(df_out)
      colnames(df_out) <- c("t_pre", "t_post", "conc_pre", "conc_post", "ss_auc")
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
