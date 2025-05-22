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

    #' Draws in a raw array with same shape as Stan variable
    #'
    #' @param name Param/quantity name of \code{x}
    #' @return Array with dimension \code{c(ndraws(x), dim(x))}
    draws_of = function(name) {
      posterior::draws_of(fit$draws(name), with_chains = FALSE)
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      S <- self$num_draws()
      x1 <- paste0("A MultistateModelStanFit with ", S, " draws.")
      msg <- paste(x1, "\n", sep = "\n")
      cat(msg)
    },

    #' @description Plot used basis functions (grid)
    plot_basisfun = function() {
      sd <- self$stan_data
      t <- rep(sd$t_grid, sd$N_sbf)
      y <- as.vector(sd$SBF_grid)
      idx <- as.factor(rep(1:sd$N_sbf, each = sd$N_grid))
      df <- data.frame(t, y, idx)
      ggplot(df, aes(.data$t, .data$y, color = .data$idx)) +
        geom_line() +
        ggtitle("Basis functions")
    },

    #' @description
    #' Full names of parameters that start with \code{log_z_}.
    log_z_pars = function() {
      nams <- names(self$draws())
      match <- grepl(nams, pattern = "log_z_")
      nams[which(match)]
    },

    #' @description Get number of draws
    num_draws = function() {
      posterior::ndraws(self$draws("lp__"))
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
      stop("Not implemented")
    }
  )
)

# Log baseline hazard distribution at times t
msmsf_log_baseline_hazard <- function(fit, t = NULL) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  sys <- fit$model$system
  if (is.null(t)) {
    t <- seq(0, sys$get_tmax(), length.out = 30)
  }
  checkmate::assert_numeric(t, min.len = 2)
  SBF <- sys$basisfun_matrix(t)
  w <- fit$draws_of("weights") # dim = c(S, H, W)
  log_w0 <- fit$draws_of("log_w0") # dim = c(S, H)
  S <- fit$num_draws()
  N <- length(t)
  bh <- matrix(0, S, N)
  for (s in seq_len(S)) {
    bh[s, ] <- sys$log_baseline_hazard(NULL, log_w0[s, 1, 1], w[s, 1, ], SBF)
  }
  list(t = t, log_h0 = bh)
}

# Spline weights
draws_df_weights <- function(fit) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  a <- as.vector(posterior::merge_chains(fit$draws("weights")))
  S <- fit$model$system$num_trans()
  W <- fit$model$system$num_weights()
  trans_idx <- rep(1:S, times = W)
  weight_idx <- rep(1:W, each = S)
  data.frame(trans_idx = trans_idx, weight_idx = weight_idx, weight = a)
}

# As draws array with single chain
draws_array_merged <- function(fit, name) {
  a <- posterior::as_draws_array(posterior::merge_chains(fit$draws(name)))
  class(a) <- "array"
  a
}
