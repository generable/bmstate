#' Create a multistate model
#'
#' @export
#' @param tm A \code{\link{TransitionMatrix}}
#' @param hazard_covs Covariates that affect the hazard. A character vector.
#' @param pk_covs Covariates that affect the PK parameters. A list with
#' elements \code{ka} \code{CL}, and \code{V2}. If \code{NULL}, a PK model
#' will not be created.
#' @param compile compile 'Stan' model?
#' @return A \code{\link{MultistateModel}} object.
create_msm <- function(tm, hazard_covs, pk_covs = NULL,
                       compile = TRUE) {
  mss <- MultistateSystem$new(tm)
  if (!is.null(pk_covs)) {
    pk <- PKModel$new(pk_covs)
  } else {
    pk <- NULL
  }
  MultistateModel$new(mss, hazard_covs, pk, compile)
}

#' Main model class
#'
#' @export
#' @field system A \code{\link{MultistateSystem}}
#' @field pk_model A \code{\link{PKModel}} or NULL
MultistateModel <- R6::R6Class("MultistateModel",

  # PRIVATE
  private = list(
    stan_model = NULL,
    hazard_covariates = NULL,
    simulate_log_hazard_multipliers = function(df_subjects, beta) {
      ts <- self$target_states()
      x <- self$covs()
      B <- length(ts)
      K <- length(x)
      checkmate::assert_matrix(beta, nrows = B, ncols = K)
      N <- nrow(df_subjects)
      S <- self$system$num_trans()
      out <- matrix(0, N, S)
      tf <- mod$system$tm()$trans_df()
      X <- df_subjects |> dplyr::select(x)
      for (s in seq_len(S)) {
        target_state <- tf$state[s]
        idx_in_beta <- which(ts == target_state)
        if (length(idx_in_beta) != 1) {
          stop("error")
        }
        beta_s <- beta[idx_in_beta, ]
        for (n in seq_len(N)) {
          out[n, s] <- sum(as.numeric(X[n, ]) * beta_s)
        }
      }
      out
    }
  ),

  # PUBLIC
  public = list(
    system = NULL,
    pk_model = NULL,

    #' @description
    #' Create model
    #'
    #' @param system A \code{\link{MultistateSystem}}
    #' @param covariates The names of the hazard covariates.
    #' @param pk_model A \code{\link{PKModel}} or NULL.
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(system, covariates, pk_model = NULL, compile = TRUE) {
      checkmate::assert_character(covariates)
      checkmate::assert_class(system, "MultistateSystem")
      if (!is.null(pk_model)) {
        checkmate::assert_class(pk_model, "PKModel")
      }
      private$hazard_covariates <- covariates
      self$pk_model <- pk_model
      self$system <- system
      if (compile) {
        self$compile()
      }
    },

    #' @description Set knot locations based on event times
    #'
    #' @param t_max Max time
    #' @param t_event Occurred event times
    #' @param num_knots Total number of knots
    set_knots = function(t_max, t_event, num_knots) {
      checkmate::assert_integerish(num_knots, len = 1, lower = 3)
      checkmate::assert_number(t_max, lower = 0)
      checkmate::assert_numeric(t_event, lower = 0, upper = t_max, min.len = 3)
      tk <- place_internal_knots(t_max, num_knots - 2, t_event)
      knots <- c(0, tk, t_max)
      self$system$set_knots(knots)
    },

    #' @description Get names of the states
    get_states = function() {
      self$system$tm()$states
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      covs <- self$covs()
      x1 <- paste0("A MultistateModel with:")
      x2 <- paste0(" - Hazard covariates: {", paste0(covs, collapse = ", "), "}")
      msg <- paste(x1, x2, "\n", sep = "\n")
      cat(msg)
      print(self$system)
      if (!is.null(self$pk_model)) {
        print(self$pk_model)
      }
    },

    #' @description Get the hazard covariates.
    covs = function() {
      private$hazard_covariates
    },

    #' @description Simulate data using the multistate model.
    #'
    #' @param N_subject number of subjects
    #' @param beta Covariate effects on each transition type.
    #' Matrix of shape \code{num_target_states} x \code{num_covs}.
    #' If \code{NULL}, a matrix of zeros is used.
    #' @param log_w0 Baseline hazard rate for all transitions
    #' @return a \code{tibble}
    simulate_data = function(N_subject = 100, beta = NULL, log_w0 = -4) {
      df_sub <- self$simulate_subjects(N_subject)
      checkmate::assert_number(log_w0)
      log_w0_vec <- rep(log_w0, self$system$num_trans())
      if (is.null(beta)) {
        L <- length(self$target_states())
        K <- length(self$covs())
        beta <- matrix(0, L, K)
      }
      self$simulate_events(df_sub, beta, log_w0_vec)
    },

    #' @description Simulate subject data.
    #'
    #' @param N_subject number of subjects
    #' @return a \code{tibble}
    simulate_subjects = function(N_subject = 100) {
      subject_id <- paste0("sim-", seq_len(N_subject))
      covs <- self$covs()
      N_covs <- length(covs)
      A <- matrix(rnorm(N_subject * N_covs), N_subject, N_covs)
      df <- data.frame(A)
      colnames(df) <- covs
      df$subject_id <- subject_id
      df$subject_idx <- seq_len(N_subject)
      as_tibble(df)
    },

    #' @description Simulate events data.
    #'
    #' @param df_subjects The subjects data frame
    #' @param beta Matrix of shape \code{num_target_states} x \code{num_covs}
    #' @param log_w0 Baseline log hazard rate, vector with length
    #' \code{num_trans}
    #' @param w_scale scale of spline weights variation
    #' @return a \code{tibble}
    simulate_events = function(df_subjects, beta, log_w0, w_scale = 0.1) {
      dt <- 1
      N <- nrow(df_subjects)
      S <- self$system$num_trans()
      L <- self$system$num_weights()
      checkmate::assert_numeric(log_w0, len = S)
      checkmate::assert_number(w_scale, lower = 0)
      w <- array(w_scale * rnorm(N * S * L), dim = c(N, S, L))
      log_w0 <- matrix(rep(log_w0, N), N, S, byrow = TRUE)
      log_m <- private$simulate_log_hazard_multipliers(df_subjects, beta)
      paths <- mod$system$simulate(w, log_w0, log_m, dt = dt)
      as_tibble(paths)
    },

    #' Get indices of states that are not source states
    #'
    #' @return integer
    target_states = function() {
      df <- self$system$tm()$states_df() |> dplyr::filter(!.data$source)
      df$state_idx
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

# Create internal knots based on event time quantiles
place_internal_knots <- function(t_max, num_knots, t_event) {
  h <- 1 / (num_knots + 1)
  knots <- stats::quantile(t_event, probs = seq(0, 1, h))
  knots[2:(length(knots) - 1)]
}
