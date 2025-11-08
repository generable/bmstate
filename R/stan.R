#' Create the main 'Stan' model
#'
#' @export
#' @param ... Arguments passed to 'CmdStanR'
#' @param filepath Deprecated.
#' @description Uses the Stan code at
#' at \code{\link{default_stan_filepath}} or a custom file if it has been
#' set using \code{options(bmstate_stan_file = ...)}.
#' @family Stan-related functions
create_stan_model <- function(filepath = NULL, ...) {
  msg <- paste0(
    "filepath argument is deprecated, set is using ",
    "options(bmstate_stan_file = FILEPATH) instead"
  )
  if (!is.null(filepath)) {
    stop(msg)
  }
  filepath <- getOption("bmstate_stan_file", default = default_stan_filepath())
  message("using stan file at ", filepath)
  cmdstanr::cmdstan_model(filepath, ...)
}

#' Get path to default Stan file
#'
#' @export
#' @family Stan-related functions
default_stan_filepath <- function() {
  fn <- "msm.stan"
  system.file(file.path("stan", fn), package = "bmstate")
}

#' Expose 'Stan' functions if they are not yet exposed
#'
#' @export
#' @param ... Passed to \code{\link{create_stan_model}}
#' @return logical telling whether they were already exposed
#' @family Stan-related functions
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
#' @inheritParams create_stan_model
#' @param method Must be one of \code{"sample"} (default),
#' \code{"pathfinder"} or \code{"optimize"}.
#' @param ... Arguments passed to the \code{sample},
#' \code{pathfinder} or \code{optimize}
#' method of the 'CmdStanR' model.
#' @return A \code{\link{MultistateModelFit}} object.
#' @family Stan-related functions
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
