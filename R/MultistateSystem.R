#' A multistate system with proportional spline hazards
#'
#' @export
MultistateSystem <- R6::R6Class("MultistateSystem",

  # PRIVATE
  private = list(
    knots = NULL,
    spline_k = 3, # degree + 1
    transmat = NULL,

    # Generate a path starting from time 0
    generate_path = function(w, log_w0, log_m, t_max, init_state) {
      # Setup
      states <- init_state
      times <- 0
      tidx <- 0
      j <- 0

      # Generate transitions
      while (times[j + 1] < t_max) {
        j <- j + 1
        t_cur <- times[j]

        # Generate next transition
        trans <- private$generate_transition(
          states[j], t_cur, t_max, w, log_w0, log_m
        )
        t_next <- trans$t

        # Update time and state
        times <- c(times, t_next)
        tidx <- c(tidx, trans$idx)
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
      is_censor <- rep(0, L)
      is_censor[L] <- 1 # last time point is always censoring time
      is_event <- rep(1, L)
      is_event[1] <- 0 # initial state is never an event
      is_event[L] <- 0 # last state is never an event (it is censoring time)
      states[L] <- states[L - 1]
      cbind(
        time = times,
        state = states,
        is_event = is_event,
        is_censor = is_censor,
        trans_idx = tidx
      )
    },

    # Generate transition given state and transition intensity functions
    generate_transition = function(state, t_init, t_max, w, log_w0, log_m) {
      tol <- 1.03
      possible <- private$transmat$possible_transitions_from(state)
      S <- private$transmat$num_trans()
      UB <- rep(1, S) # TODO: better upper bound
      possible <- possible[which(UB[possible] > 1e-9)] # so rare that not possible
      if (length(possible) == 0) {
        # Absorbing state, no transitions possible
        return(list(t = t_max, new_state = 0, idx = 0))
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
        draw <- private$draw_time_cnhp(
          w[j, ], log_w0[j], log_m[j], t_min_found, t_init, UB[j]
        )
        if (draw$t < t_min_found) {
          trans_idx <- possible[j]
          t_min_found <- draw$t
        }
      }
      if (trans_idx > 0) {
        new_state <- private$transmat$target_state(trans_idx)
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
        t = t_min_found, new_state = new_state, idx = trans_idx
      )
    },

    # Draw an event time from Censored Non-Homogeneous Poisson
    # process using the thinning algorithm
    draw_time_cnhp = function(w, log_w0, log_m, t_max, t_init, lambda_ub) {
      checkmate::assert_number(lambda_ub, lower = 0)
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
      # Return
      list(t = t, n_tries = n_tries, censor = FALSE)
    }
  ),

  # PUBLIC
  public = list(

    #' @description
    #' Create model
    #'
    #' @param transmat A \code{\link{TransitionMatrix}}.
    initialize = function(transmat) {
      checkmate::assert_class(transmat, "TransitionMatrix")
      private$transmat <- transmat
    },

    #' @description Get number of states
    #'
    #' @return integer
    num_states = function() {
      self$tm()$num_states()
    },

    #' @description Get number of transitions
    #'
    #' @return integer
    num_trans = function() {
      self$tm()$num_trans()
    },

    #' Get number of spline weight parameters
    #' @return an integer (\code{num_knots} - 2) + \code{self$spline_k}
    num_weights = function() {
      if (self$knots_not_set()) {
        stop("knots have not been set, use $set_knots()")
      }
      self$num_knots() - 2 + private$spline_k
    },

    #' @description Get number of spline knots
    #'
    #' @return integer, zero if knots have not been set
    num_knots = function() {
      if (self$knots_not_set()) {
        return(0)
      }
      length(private$knots)
    },

    #' @description Get transition matrix
    tm = function() {
      private$transmat
    },

    #' @description Have knots not been set yet?
    #'
    #' @return a logical
    knots_not_set = function() {
      is.null(private$knots)
    },

    #' @description Set spline knot locations
    #'
    #' @param locations A numeric vector
    set_knots = function(locations) {
      checkmate::assert_numeric(locations, min.len = 2)
      private$knots <- locations
    },

    #' @description Get knot locations
    #'
    #' @return a numeric vector
    get_knots = function() {
      if (self$knots_not_set()) {
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
      s <- self$tm()$states
      x1 <- paste0("A MultistateSystem with:")
      x2 <- paste0(" - States: {", paste0(s, collapse = ", "), "}")
      x3 <- paste0(" - Number of spline knots: ", self$num_knots())
      msg <- paste(x1, x2, x3, "\n", sep = "\n")
      cat(msg)
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
      knots <- self$get_knots()
      L <- length(knots)
      BK <- knots[c(1, L)]
      knots <- knots[2:(L - 1)]
      SBF <- bspline_basis(t, private$spline_k, knots, BK) # shape (N, L+1)
      if (is.null(w)) {
        w <- rep(0, self$num_weights())
      }
      as.numeric(SBF %*% w + log_w0)
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
    #' @return a data frame
    simulate = function(w, log_w0, log_m, init_state = 1) {
      checkmate::assert_array(w, d = 3)
      n_paths <- dim(w)[1]
      checkmate::assert_true(dim(w)[3] == self$num_weights())
      checkmate::assert_matrix(log_w0, nrows = n_paths)
      checkmate::assert_matrix(log_m, nrows = n_paths)
      S <- self$num_states()
      checkmate::assert_integerish(init_state, len = 1, lower = 1, upper = S)
      pb <- progress::progress_bar$new(total = n_paths)
      message("Generating ", n_paths, " paths")

      # Could be done in parallel and some things be precomputed
      out <- NULL
      t_max <- self$get_tmax()
      for (j in seq_len(n_paths)) {
        pb$tick()
        p <- private$generate_path(
          w[j, , ], log_w0[j, ], log_m[j, ], t_max, init_state
        )
        p <- cbind(p, rep(j, nrow(p)))
        out <- rbind(out, p)
      }
      df <- data.frame(out)
      colnames(df)[ncol(df)] <- "path_id"
      df
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
