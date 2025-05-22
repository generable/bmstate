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
      posterior::draws_of(self$draws(name), with_chains = FALSE)
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

    #' Plot baseline hazard distribution
    #'
    #' @param t times where to evaluate the baseline hazards
    #' @param ci_alpha width of credible interval
    plot_h0 = function(t = NULL, ci_alpha = 0.95) {
      df <- self$h0_dist(t, ci_alpha)
      legend <- self$model$system$tm()$trans_df()
      df <- df |> dplyr::left_join(legend, by = "trans_idx")
      df$trans_type <- as.factor(df$trans_type)
      capt <- paste0("Median and ", 100 * ci_alpha, "% CrI")
      ggplot(df, aes(
        x = .data$time, y = .data$median, ymin = .data$lower,
        ymax = .data$upper, color = .data$trans_type, fill = .data$trans_type
      )) +
        geom_ribbon(alpha = 0.7) +
        geom_line() +
        facet_wrap(. ~ .data$trans_char) +
        ylab("Log baseline hazard") +
        theme(legend.position = "none") +
        labs(caption = capt)
    },

    #' Baseline hazard distribution
    #'
    #' @param t times where to evaluate the baseline hazards
    #' @param ci_alpha width of credible interval
    h0_dist = function(t = NULL, ci_alpha = 0.95) {
      checkmate::assert_number(ci_alpha, lower = 0, upper = 1)
      LB <- (1 - ci_alpha) / 2
      UB <- 1 - LB
      df <- msmsf_log_baseline_hazard(self, t)
      df |>
        dplyr::group_by(.data$time, .data$trans_idx) |>
        dplyr::summarize(
          median = stats::median(.data$log_h0),
          upper = stats::quantile(.data$log_h0, UB),
          lower = stats::quantile(.data$log_h0, LB)
        ) |>
        dplyr::ungroup()
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

msmsf_log_m_per_subject <- function(fit) {
  log_m <- fit$draws_of("log_C_haz")
  idx_sub <- fit$stan_data$idx_sub
  N_sub <- fit$stan_data$N_sub
  first_indices <- sapply(seq_len(N_sub), function(x) which(idx_sub == x)[1])
  log_m[, first_indices, ]
}

msmsf_pathgen <- function(fit, init_state = 1) {
  checkmate::assert_class(fit, "MultistateModelStanFit")
  sys <- fit$model$system
  S <- fit$num_draws()
  N <- fit$stan_data$N_sub
  w <- fit$draws_of("weights")
  log_w0 <- fit$draws_of("log_w0")
  w_rep <- abind::abind(replicate(N, w, simplify = FALSE), along = 1)
  log_w0_rep <- abind::abind(replicate(N, log_w0, simplify = FALSE), along = 1)
  log_m <- msmsf_log_m_per_subject(fit)
  H <- dim(log_m)[3]
  log_m_reshaped <- matrix(aperm(log_m, c(1, 2, 3)), nrow = N * S, ncol = H)
  sim <- sys$simulate(w_rep, log_w0_rep, log_m_reshaped, init_state)
  sim
}
