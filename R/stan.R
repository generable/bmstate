# Checks to be made that make sure that the model settings
# make sense for given data
prefit_checks <- function(model, data) {
  tmax <- model$get_tmax()
  pd <- data$paths
  lens <- pd$as_transitions() |>
    dplyr::mutate(time_len = .data$time - .data$time_prev) |>
    dplyr::pull(.data$time_len)
  min_len <- min(lens)
  delta_grid <- tmax / model$n_grid
  if (min_len < delta_grid) {
    msg <- paste0(
      "Shortest time interval (", min_len,
      ") is smaller than delta_grid (", delta_grid,
      "). Consider increasing n_grid or decreasing t_max of the model."
    )
    message(msg)
  }
  max_time <- max(pd$get_path_df()$time)
  if (tmax < max_time) {
    msg <- paste0(
      "Model t_max (", tmax,
      ") is smaller than max observed time in data (",
      max_time, "). Increase t_max."
    )
    stop(msg)
  }
  if (tmax > 2 * max_time) {
    warning(
      "Set model t_max is more than twice as large as largest observed time",
      " in the data. Are you sure you want to do this? If not, set smaller",
      " t_max."
    )
  }
  TRUE
}

#' Create the main 'Stan' model
#'
#' @export
#' @param ... Arguments passed to 'CmdStanR'
#' @param filepath File path to the \code{msm.stan} file.
create_stan_model <- function(filepath = NULL, ...) {
  fn <- "msm.stan"
  if (is.null(filepath)) {
    filepath <- system.file(file.path("stan", fn), package = "bmstate")
  }
  cmdstanr::cmdstan_model(filepath, ...)
}

#' Expose 'Stan' functions if they are not yet exposed
#'
#' @export
#' @param ... Passed to \code{\link{create_stan_model}}
#' @return logical telling whether they were already exposed
ensure_exposed_stan_functions <- function(...) {
  if (exists("STAN_dummy_function")) {
    return(TRUE)
  } else {
    message("Recompiling Stan model")
    mod <- create_stan_model(force_recompile = TRUE, ...)
    mod$expose_functions(global = TRUE)
  }
  FALSE
}

#' Fit a model using 'Stan'
#'
#' @description
#' \emph{NOTE:} This function has a side effect of setting covariate
#' normalizers and prior assumed mean baseline hazard based on data.
#'
#' @export
#' @param model A \code{\link{MultistateModel}} object.
#' @param data A \code{\link{JointData}} or \code{\link{PathData}} object.
#' @param return_stanfit Return also the raw 'Stan' fit object?
#' @param set_auc_normalizers Set AUC normalization based on SS doses.
#' @param filepath Passed to \code{\link{create_stan_model}}.
#' @param method Must be one of \code{"sample"} (default),
#' \code{"pathfinder"} or \code{"optimize"}.
#' @param ... Arguments passed to the \code{sample},
#' \code{pathfinder} or \code{optimize}
#' method of the 'CmdStanR' model.
#' @return A \code{\link{MultistateModelFit}} object.
fit_stan <- function(model, data,
                     set_auc_normalizers = TRUE,
                     filepath = NULL,
                     return_stanfit = FALSE,
                     method = "sample", ...) {
  checkmate::assert_class(model, "MultistateModel")
  data <- pd_to_jointdata(data)
  checkmate::assert_class(data, "JointData")
  checkmate::assert_character(method, len = 1, min.chars = 1)
  checkmate::assert_choice(method, c("sample", "pathfinder", "optimize"))

  checkmate::assert_logical(return_stanfit, len = 1)
  checkmate::assert_logical(set_auc_normalizers, len = 1)
  prefit_checks(model, data)

  # Get Stan model object
  stan_model <- create_stan_model(filepath = filepath)

  # Set normalizing locations and scales (side effect)
  model$set_normalizers(data)
  if (!is.null(data$dosing) && set_auc_normalizers) {
    mu_CL <- exp(-2)
    aaa <- data$dosing$dose_ss / mu_CL
    loc <- mean(aaa)
    sca <- stats::sd(aaa)
    model$set_auc_normalizers(loc, sca)
  }

  # Set prior mean baseline hazard rates (side effect)
  df_ttype <- average_haz_per_ttype(data$paths) |> dplyr::arrange(.data$trans_idx)
  model$set_prior_mean_h0(exp(df_ttype$log_h0_avg))

  # Create Stan input list
  sd <- create_stan_data(model, data)

  # Call 'Stan'
  if (method == "pathfinder") {
    stan_fit <- stan_model$pathfinder(data = sd, ...)
  } else if (method == "optimize") {
    stan_fit <- stan_model$optimize(data = sd, ...)
  } else {
    stan_fit <- stan_model$sample(data = sd, ...)
  }

  # Return
  pars <- c(
    "weights", "log_w0", "beta_ka", "beta_V2", "beta_CL", "beta_oth",
    "beta_auc", "sigma_pk", "log_z_pk", "log_mu_pk", "log_sig_pk", "lp__"
  )
  draws <- create_rv_list(stan_fit, pars)
  diag <- NULL
  if (method == "sample") {
    diag <- stan_fit$diagnostic_summary()
  }
  sfun <- function(x) x$summary()
  quiet_sumr <- purrr::quietly(sfun)
  smr <- quiet_sumr(stan_fit)
  info <- list(
    diag = diag,
    runset = stan_fit$runset,
    summary = smr$result
  )
  fit <- MultistateModelFit$new(data, sd, model, draws, info)
  if (!return_stanfit) {
    return(fit)
  }
  list(
    fit = fit,
    stanfit = stan_fit
  )
}
