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

    #' @description Evaluate log instant hazard
    #'
    #' @param t Time point(s)
    #' @param w Spline basis function weights (vector)
    #' @param log_w0 Intercept (log)
    #' @param log_m Hazard multiplier (log)
    log_inst_hazard = function(t, w, log_w0, log_m) {
      log_m + self$log_baseline_hazard(t, log_w0, w)
    },

    #' Generate paths
    #'
    #' @param w An array of shape \code{n_paths} x \code{n_trans} x
    #' \code{n_weights}
    #' @param log_w0 An array of shape \code{n_paths} x \code{n_trans}
    #' @param log_m An array of shape \code{n_paths} x \code{n_trans}
    #' @param init_state Integer index of starting state.
    #' @param discretize Discretize times to one time unit?
    #' @return a data frame
    simulate = function(w, log_w0, log_m, init_state = 1, discretize = FALSE) {
      checkmate::assert_array(w, d = 3)
      n_paths <- dim(w)[1]
      checkmate::assert_matrix(log_w0, nrows = n_paths)
      checkmate::assert_matrix(log_m, nrows = n_paths)
      checkmate::assert_logical(discretize, len = 1)
      S <- self$transmat$num_states()
      checkmate::assert_integerish(init_state, len = 1, lower = 1, upper = S)
      pb <- progress::progress_bar$new(total = n_paths)
      message("Generating ", n_paths, " paths")

      # Could be done in parallel
      out <- NULL
      for (j in seq_len(n_paths)) {
        pb$tick()
        p <- self$generate_path(
          w[j, , ], log_w0[j, ], log_m[j, ], init_state, discretize
        )
        p <- cbind(p, rep(j, nrow(p)))
        out <- rbind(out, p)
      }
      out
    },

    #' Generate a path starting from time 0
    #'
    #' @param w An array of shape \code{n_trans} x \code{n_weights}
    #' @param log_w0 A vector of length \code{n_trans}
    #' @param log_m A vector of length \code{n_trans}
    #' @param init_state Integer index of starting state.
    #' @param discretize Discretize times to one time unit?
    #' @return an array with number of rows equal to the path length
    generate_path = function(w, log_w0, log_m,
                             init_state = 1, discretize = FALSE) {
      checkmate::assert_logical(discretize, len = 1)
      checkmate::assert_integerish(init_state, len = 1)
      checkmate::assert_array(w, d = 2)
      checkmate::assert_numeric(log_w0)
      checkmate::assert_numeric(log_m)

      t_max <- self$get_tmax()
      S <- self$transmat$num_states()
      checkmate::assert_integerish(init_state, len = 1, lower = 1, upper = S)
      dt <- 1

      # Setup
      states <- init_state
      times <- 0
      j <- 0

      # Generate transitions
      while (times[j + 1] < t_max) {
        j <- j + 1
        t_cur <- times[j]

        # Generate next transition
        trans <- self$generate_transition(
          states[j], t_cur, t_max, w, log_w0, log_m
        )
        t_next <- trans$t

        # Discretize transition time to end of day
        if (discretize) {
          t_next <- dt * ceiling(t_next / dt)
        }

        # Update time and state
        times <- c(times, t_next)
        states <- c(states, trans$new_state)
        if (j >= 100000) {
          msg <- paste0(
            "\nt = ",
            round(times[j], 2), " -> ", t_max, ": at state ", states[j], "\n ",
            "Generating a path with over ", j, " transitions, something is wrong\n"
          )
          stop(msg)
        }
      }
      L <- length(times)
      is_event <- rep(1, L)
      is_event[1] <- 0 # initial state is never an event
      is_event[L] <- 0
      states[L] <- states[L - 1]
      cbind(time = times, state = states, is_event = is_event)
    },

    #' Generate transition given state and transition intensity functions
    #'
    #' @param state Integer index of current state.
    #' @param t_init Current time
    #' @param t_max max time
    #' @param w An array of shape \code{n_trans} x \code{n_weights}
    #' @param log_w0 A vector of length \code{n_trans}
    #' @param log_m A vector of length \code{n_trans}
    generate_transition = function(state, t_init, t_max, w, log_w0, log_m) {
      checkmate::assert_integerish(state, len = 1)
      checkmate::assert_number(t_init)
      checkmate::assert_number(t_max, lower = t_init)
      checkmate::assert_array(w, d = 2)
      checkmate::assert_numeric(log_w0)
      checkmate::assert_numeric(log_m)
      tol <- 1.03
      possible <- self$transmat$possible_transitions_from(state)
      S <- self$transmat$num_trans()
      UB <- rep(1, S) # TODO: better upper bound
      possible <- possible[which(UB[possible] > 1e-9)] # so rare that not possible
      if (length(possible) == 0) {
        # Absorbing state, no transitions possible
        return(list(t = t_max, new_state = 0))
      }
      UB <- UB[possible]
      w <- w[possible, , drop = FALSE]
      log_w0 <- log_w0[possible]
      log_m <- log_m[possible]

      # Draw event with highest upper bound first because it is most likely
      # to occur soon and be the minimum
      draw_order <- sort(UB, decreasing = TRUE, index.return = 1)$ix
      t_min_found <- t_max
      trans_idx <- 0

      # Loop through possible transitions
      for (j in draw_order) {
        draw <- self$draw_time_cnhp(
          w[j, ], log_w0[j], log_m[j], t_min_found, t_init, UB[j]
        )
        if (draw$t < t_min_found) {
          trans_idx <- possible[j]
          t_min_found <- draw$t
        }
      }
      if (trans_idx > 0) {
        new_state <- self$transmat$target_state(trans_idx)
      } else {
        new_state <- 0
      }

      # Safeguard against infinite loop
      dt_min <- 1e-9
      if (t_min_found - t_init < dt_min) {
        t_min_found <- t_init + dt_min
      }

      # Return
      list(
        t = t_min_found, new_state = new_state
      )
    },

    #' @description Draw an event time from Censored Non-Homogeneous Poisson
    #' process using the thinning algorithm
    #'
    #' @param w Spline basis function weights (vector)
    #' @param log_w0 Intercept (log)
    #' @param log_m Hazard multiplier (log)
    #' @param t_max Max time
    #' @param t_init Initial time
    #' @param lambda_ub Upper bound of hazard
    #' @return A number. If drawn event time is going to be larger than t_max,
    #' then t_max is returned
    draw_time_cnhp = function(w, log_w0, log_m, t_max, t_init, lambda_ub) {
      n_tries <- 0
      t <- t_init
      accepted <- FALSE
      while (!accepted) {
        n_tries <- n_tries + 1
        u1 <- stats::runif(1)
        t <- t - 1 / lambda_ub * log(u1)
        if (t > t_max) {
          return(list(t = t_max, n_tries = n_tries, censor = TRUE))
        }
        lambda_t <- exp(self$log_inst_hazard(t, w, log_w0, log_m))
        if (lambda_t > lambda_ub) {
          msg <- "lambda(t) was larger than its supposed upper bound"
          msg <- paste0(
            msg, "\nt = ", t, ", t_init = ", t_init,
            ", \nlambda(t) = ", lambda_t, ", UB = ", lambda_ub
          )
          stop(msg)
        }
        p <- lambda_t / lambda_ub
        u2 <- stats::runif(1)
        if (u2 <= p) {
          accepted <- TRUE
        }
      }
      return(list(t = t, n_tries = n_tries, censor = FALSE))
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
