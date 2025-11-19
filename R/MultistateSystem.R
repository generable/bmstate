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
    generate_path = function(w, log_w0, log_m, t_start, t_max, init_state, min_t_step) {
      if (is.vector(w)) {
        w <- matrix(w, nrow = 1)
      }
      checkmate::assert_array(w, d = 2)

      # Setup
      states <- init_state
      times <- t_start
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
        if (t_next - t_cur < min_t_step) {
          t_next <- t_cur + min_t_step
        }

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
      tol <- 1.05
      possible <- private$transmat$possible_transitions_from(state)
      J <- private$transmat$num_trans()
      UB <- rep(0, J)
      for (j in seq_len(J)) {
        UB[j] <- tol * self$max_inst_hazard(t_init, t_max, w[j, ], log_w0[j], log_m[j])
      }
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

    #' @description Get number of spline weight parameters
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
      names(locations) <- NULL
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

    #' @description Print the object
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
    #' @param t Output time points. Not used if \code{SBF} is given.
    #' @param log_w0 Intercept (log_scale).
    #' @param w Spline weights (log scale). If \code{NULL}, will be set to
    #' a vector of zeros, meaning that the log hazard is constant at
    #' \code{log_w0}
    #' @param SBF Precomputed basisfunction matrix at \code{t}.
    #' @return a vector with same length as \code{t}
    log_baseline_hazard = function(t, log_w0, w = NULL, SBF = NULL) {
      if (is.null(SBF)) {
        SBF <- self$basisfun_matrix(t)
      }
      if (is.null(w)) {
        w <- rep(0, self$num_weights())
      }
      as.numeric(SBF %*% w + log_w0)
    },

    #' @description Evaluate basis function matrix for baseline hazard
    #'
    #' @param t Evaluation time points.
    #' @return Matrix with shape \code{c(N, L+1)} where \code{L} is number of
    #' knots
    basisfun_matrix = function(t) {
      knots <- self$get_knots()
      L <- length(knots)
      BK <- knots[c(1, L)]
      knots <- knots[2:(L - 1)]
      bspline_basis(t, private$spline_k, knots, BK)
    },

    #' @description Evaluate log instant hazard
    #'
    #' @param t Time point(s). Not used if \code{SBF} is given.
    #' @param w Spline basis function weights (vector)
    #' @param log_w0 Intercept (log)
    #' @param log_m Hazard multiplier (log)
    #' @param SBF Pre-computed basis function matrix at \code{t}.
    log_inst_hazard = function(t, w, log_w0, log_m, SBF = NULL) {
      log_m + self$log_baseline_hazard(t, log_w0, w, SBF)
    },

    #' @description Evaluate transition intensity (generator) matrix at time t
    #'
    #' @param t A number (time point)
    #' @param w An array of shape \code{n_trans} x \code{n_weights}
    #' @param log_w0 A vector of length \code{n_trans}
    #' @param log_m A vector of length \code{n_trans}
    #' @return a matrix with shape \code{n_states} x \code{n_states}
    intensity_matrix = function(t, log_w0, w = NULL, log_m = NULL) {
      if (self$has_self_loops()) {
        stop("System is not a standard continuous-time Markov multistate model")
      }
      tm <- self$tm()
      mat <- tm$as_transition_index_matrix()
      H <- self$num_trans()
      checkmate::assert_numeric(log_w0, len = H)
      S <- self$num_states()
      if (is.null(w)) {
        W <- self$num_weights()
        w <- matrix(0, H, W)
      }
      if (is.null(log_m)) {
        log_m <- rep(0, H)
      }
      ti <- find_row_and_col_of_positive_vals(mat)
      for (idx in seq_len(H)) {
        log_h <- self$log_inst_hazard(t, w[idx, ], log_w0[idx], log_m[idx])
        mat[ti[idx, 1], ti[idx, 2]] <- exp(log_h)
      }
      for (s in seq_len(S)) {
        mat[s, s] <- -sum(mat[s, ])
      }
      mat
    },

    #' @description Does the system have self loops?
    #' @return Boolean value
    has_self_loops = function() {
      sum(diag(self$tm()$matrix)) > 0
    },

    #' @description Max instant hazard on interval (t1, t2)
    #'
    #' @param t1 Start time point
    #' @param t2 End time point
    #' @param w Spline basis function weights (vector)
    #' @param log_w0 Intercept (log)
    #' @param log_m Hazard multiplier (log)
    max_inst_hazard = function(t1, t2, w, log_w0, log_m) {
      ttt <- seq(t1, t2, length.out = 100)
      log_haz <- log_m + self$log_baseline_hazard(ttt, log_w0, w)
      return(exp(max(log_haz)))
    },

    #' @description Generate paths
    #'
    #' @param w An array of shape \code{n_draws} x \code{n_trans} x
    #' \code{n_weights}
    #' @param log_w0 An array of shape \code{n_draws} x \code{n_trans}
    #' @param log_m An array of shape \code{n_draws} x \code{n_trans}
    #' @param init_state Index of starting state. A single value or a vector
    #' with length equal to \code{n_draws}.
    #' @param t_start Start time.
    #' @param t_max Max time. If \code{NULL}, the max
    #' time of the model is used.
    #' @param n_rep Number of repetitions to do for each draw.
    #' @param min_t_step Minimal time step.
    #' @return A data frame with \code{n_draws} x \code{n_rep} paths.
    simulate = function(w, log_w0, log_m, init_state = 1, t_start = 0,
                        t_max = NULL, n_rep = 1, min_t_step = 1e-6) {
      checkmate::assert_array(w, d = 3)
      n_draws <- dim(w)[1]
      checkmate::assert_true(dim(w)[3] == self$num_weights())
      checkmate::assert_matrix(log_w0, nrows = n_draws)
      checkmate::assert_matrix(log_m, nrows = n_draws)
      checkmate::assert_number(min_t_step, lower = 0)
      S <- self$num_states()
      checkmate::assert_integerish(init_state, lower = 1, upper = S)
      if (length(init_state) == 1) {
        init_state <- rep(init_state, n_draws)
      } else {
        stopifnot(length(init_state) == n_draws)
      }
      checkmate::assert_number(t_start, lower = 0)
      checkmate::assert_integerish(n_rep, len = 1, lower = 1)
      n_paths <- n_draws * n_rep
      pb <- progress::progress_bar$new(total = n_paths)

      # Set max time
      out <- NULL
      if (is.null(t_max)) {
        t_max <- self$get_tmax()
      }
      checkmate::assert_number(t_max, lower = 0)

      # Should not be done in parallel as such because can mess order in link df
      cnt <- 0
      message("Generating ", n_paths, " paths")
      for (k in seq_len(n_rep)) {
        for (j in seq_len(n_draws)) {
          cnt <- cnt + 1
          pb$tick()
          p <- private$generate_path(
            w[j, , ], log_w0[j, ], log_m[j, ], t_start, t_max, init_state[j],
            min_t_step
          )
          p <- cbind(p, rep(cnt, nrow(p)))
          out <- rbind(out, p)
        }
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
    degree = k - 1, knots = knots, intercept = TRUE, Boundary.knots = BK
  )
}

#' Solve the Kolmogorov forward equation of a Markovian multistate system
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param t A time grid
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_w0 A vector of length \code{n_trans}
#' @param log_m A vector of length \code{n_trans}
#' @param ... Arguments passed to \code{deSolve::ode()}.
#' @return Value given by \code{deSolve::ode()}.
solve_time_evolution <- function(system, t, log_w0, w = NULL, log_m = NULL,
                                 ...) {
  checkmate::assert_class(system, "MultistateSystem")
  S <- system$num_states()
  H <- system$num_trans()
  W <- system$num_weights()
  if (is.null(w)) {
    w <- matrix(0, H, W)
  }
  if (is.null(log_m)) {
    log_m <- rep(0, H)
  }
  checkmate::assert_numeric(t, min.len = 1)
  checkmate::assert_matrix(w, ncols = W, nrows = H)
  checkmate::assert_numeric(log_w0, len = H)
  checkmate::assert_numeric(log_m, len = H)

  odefun <- function(time, y, parms) {
    P <- matrix(y, S, S)
    Lambda <- system$intensity_matrix(time, log_w0, w, log_m)
    dydt <- as.vector(P %*% Lambda)
    list(dydt)
  }
  P0 <- diag(1, S, S)
  y0 <- as.vector(P0)
  deSolve::ode(y0, t, odefun, NULL, method = "ode45", ...)
}

#' Solve the transition probability matrix for each given output time
#'
#' @export
#' @param system A \code{\link{MultistateSystem}}
#' @param log_w0 A vector of length \code{n_trans}
#' @param t_start Initial time
#' @param t_out End times (vector)
#' @param w An array of shape \code{n_trans} x \code{n_weights}
#' @param log_m A vector of length \code{n_trans}
#' @param ... Arguments passed to \code{deSolve::ode()}.
#' @return An array \code{P} where \code{P[k,i,j]} is the probability that
#' the system will be in state \code{j} at the \code{k}th time point of
#' \code{t_out} given that it is in state \code{i} at time \code{t_start}.
solve_trans_prob_matrix <- function(system, t_out, log_w0, w = NULL,
                                    log_m = NULL, t_start = 0,
                                    ...) {
  checkmate::assert_class(system, "MultistateSystem")
  checkmate::assert_numeric(t_out)
  K <- length(t_out)
  checkmate::assert_number(t_start, lower = 0)
  H <- system$num_trans()
  W <- system$num_weights()
  S <- system$num_states()
  if (is.null(w)) {
    w <- matrix(0, H, W)
  }
  if (is.null(log_m)) {
    log_m <- rep(0, H)
  }
  kfe <- solve_time_evolution(system, t_out, log_w0, w, log_m, ...)
  P <- array(0, dim = c(K, S, S))
  for (k in seq_len(K)) {
    P[k, , ] <- matrix(kfe[k, 2:ncol(kfe)], S, S)
  }
  P
}
