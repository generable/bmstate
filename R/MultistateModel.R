#' Create a multistate model
#'
#' @export
#' @param tm A \code{\link{TransitionMatrix}}. See \code{\link{transmat}}
#' for how to create common transition matrices.
#' @param hazard_covs Covariates that affect the hazard. A character vector.
#' The name \code{"dose_amt"} is special if simulating data using the model.
#' @param categ_covs Names of covariates that are binary.
#' This only has an effect when simulating data. When fitting a model, all
#' covariates are treated as continuous, so you should use a binary encoding
#' for categories if there is more than two.
#' @param pk_covs Covariates that affect the PK parameters. A list with
#' elements \code{ka} \code{CL}, and \code{V2}. If \code{NULL}, a PK model
#' will not be created.
#' @param ... Arguments passed to \code{\link{MultistateModel}} init
#' @return A \code{\link{MultistateModel}} object.
create_msm <- function(tm, hazard_covs = NULL, pk_covs = NULL,
                       categ_covs = NULL, ...) {
  mss <- MultistateSystem$new(tm)
  if (!is.null(pk_covs)) {
    pk <- PKModel$new(pk_covs)
  } else {
    pk <- NULL
  }
  MultistateModel$new(mss, hazard_covs, pk,
    categorical = categ_covs,
    ...
  )
}

#' Main model class
#'
#' @description Class that represents a multistate model.
#'
#' @export
#' @field system A \code{\link{MultistateSystem}}
#' @field pk_model A \code{\link{PKModel}} or NULL.
#' @field prior_only Should the model ignore likelihood?
#' @field pk_only Should the model ignore the entire
#' hazard model part?
MultistateModel <- R6::R6Class("MultistateModel",

  # PRIVATE
  private = list(
    prior_mean_h0 = NULL,
    hazard_covariates = NULL,
    categorical = NULL,
    normalizer_locations = NULL,
    normalizer_scales = NULL,
    auc_normalizer_loc = 300,
    auc_normalizer_scale = 100,
    n_grid = NULL,
    simulate_log_hazard_multipliers = function(df_subjects, beta) {
      ts <- self$target_states()
      x <- self$covs()
      auc_norm <- self$get_auc_normalizers()
      B <- length(ts)
      K <- length(x)
      checkmate::assert_matrix(beta, nrows = B, ncols = K)
      N <- nrow(df_subjects)
      S <- self$system$num_trans()
      out <- matrix(0, N, S)
      tf <- self$system$tm()$trans_df()
      X <- df_subjects |> dplyr::select(tidyselect::all_of(x))
      X <- as.matrix(X)
      X_norm <- normalize_columns(X)
      if ("ss_auc" %in% x) {
        idx <- which(x == "ss_auc")
        x_norm <- (X[, idx] - auc_norm$loc) / auc_norm$scale
        check_normalized_covariate(x_norm, "ss_auc")
        X_norm[, idx] <- x_norm
      }
      for (s in seq_len(S)) {
        target_state <- tf$state[s]
        idx_in_beta <- which(ts == target_state)
        if (length(idx_in_beta) != 1) {
          stop("error")
        }
        beta_s <- beta[idx_in_beta, ]
        for (n in seq_len(N)) {
          out[n, s] <- sum(as.numeric(X_norm[n, ]) * beta_s)
        }
      }
      out
    },

    # Simulate PK data
    simulate_pk_data = function(df_subjects, beta_pk = NULL) {
      pk_dat <- NULL
      if (self$has_pk()) {
        beta_pk <- self$pk_model$format_params(beta_pk)
        pk_dat <- self$pk_model$simulate_data(df_subjects, beta_pk)
      }
      pk_dat
    },

    # Simulate events data
    simulate_events = function(df_subjects, beta_haz, log_w0, w) {
      dt <- 1
      N <- nrow(df_subjects)
      S <- self$system$num_trans()
      L <- self$system$num_weights()
      checkmate::assert_numeric(log_w0, len = S)
      w_all <- array(0, dim = c(N, S, L))
      if (!is.null(w)) {
        for (n in seq_len(N)) {
          w_all[n, , ] <- w
        }
      }

      log_w0 <- matrix(rep(log_w0, N), N, S, byrow = TRUE)
      log_m <- private$simulate_log_hazard_multipliers(df_subjects, beta_haz)
      paths <- self$system$simulate(w_all, log_w0, log_m, min_t_step = 0.1)
      as_tibble(paths)
    }
  ),

  # PUBLIC
  public = list(
    system = NULL,
    pk_model = NULL,
    prior_only = FALSE,
    pk_only = FALSE,

    #' @description Get normalization constants for each variable
    #' @return list
    get_normalizers = function() {
      list(
        locations = private$normalizer_locations,
        scales = private$normalizer_scales
      )
    },

    #' @description Get number of grid points used for integration.
    #' @return An integer
    get_n_grid = function() {
      private$n_grid
    },

    #' @description Set normalization constant for each variable (side effect)
    #'
    #' @param data A \code{\link{JointData}} object
    set_normalizers = function(data) {
      checkmate::assert_class(data, "JointData")
      df_sub <- data$paths$subject_df
      num_cols <- which(sapply(df_sub, is.numeric))
      private$normalizer_locations <- lapply(df_sub[, num_cols], mean)
      private$normalizer_scales <- lapply(df_sub[, num_cols], stats::sd)
      invisible(NULL)
    },

    #' @description Get normalization constants for AUC (PK)
    #' @return list
    get_auc_normalizers = function() {
      list(
        loc = private$auc_normalizer_loc,
        scale = private$auc_normalizer_scale
      )
    },

    #' @description Set normalization constants for AUC (side effect)
    #'
    #' @param loc Location
    #' @param scale Scale
    set_auc_normalizers = function(loc = 0, scale = 1) {
      checkmate::assert_numeric(loc, lower = 0, len = 1)
      checkmate::assert_numeric(scale, lower = 0, len = 1)
      message(
        "setting auc normalizers to loc = ",
        round(loc, 5), ", scale = ", round(scale, 5)
      )
      private$auc_normalizer_loc <- loc
      private$auc_normalizer_scale <- scale
      invisible(NULL)
    },

    #' @description Get assumed prior mean baseline hazard rates.
    #' @return Numeric vector with length equal to number of transitions
    get_prior_mean_h0 = function() {
      v <- private$prior_mean_h0
      if (is.null(v)) {
        stop("prior mean h0 has not been set")
      }
      v
    },

    #' @description Set assumed prior mean baseline hazard rates (side
    #' effect).
    #' @param mean_h0 Numeric vector with length equal to number of transitions
    set_prior_mean_h0 = function(mean_h0) {
      N_trans <- self$system$tm()$num_trans()
      checkmate::assert_numeric(mean_h0, len = N_trans, lower = 0)
      private$prior_mean_h0 <- mean_h0
      invisible(NULL)
    },

    #' @description Set assumed prior mean baseline hazard rates (side
    #' effect) based on average hazards in data.
    #' @param data A \code{\link{JointData}} or \code{\link{PathData}} object.
    set_prior_mean_h0_data = function(data) {
      if (inherits(data, "PathData")) {
        data <- JointData$new(data, NULL)
      }
      checkmate::assert_class(data, "JointData")
      df_ttype <- average_haz_per_ttype(data$paths) |> dplyr::arrange(.data$trans_idx)
      self$set_prior_mean_h0(exp(df_ttype$log_h0_avg))
    },

    #' @description Create model
    #'
    #' @param system A \code{\link{MultistateSystem}}
    #' @param covariates The names of the hazard covariates (excluding possible
    #' exposure estimated from PK model). Do not use reserved names
    #' \code{ss_auc} or \code{dose}.
    #' @param pk_model A \code{\link{PKModel}} or NULL.
    #' @param t_max Max time.
    #' @param num_knots Total number of spline knots.
    #' @param categorical Names of covariates that are binary.
    #' This only has an effect when simulating data. When fitting a model, all
    #' covariates are treated as continuous, so you should use a binary encoding
    #' for categories if there is more than two.
    #' @param n_grid Number of time discretization points for integrating
    #' @param prior_only Should the model ignore likelihood?
    #' @param pk_only Should the model ignore the entire
    #' hazard model part?
    #' hazards.
    initialize = function(system, covariates = NULL, pk_model = NULL,
                          t_max = 1000, num_knots = 5, categorical = NULL,
                          n_grid = 1000, prior_only = FALSE, pk_only = FALSE) {
      checkmate::assert_character(covariates, null.ok = TRUE)
      checkmate::assert_character(categorical, null.ok = TRUE)
      if (!all(categorical %in% covariates)) {
        stop("all categorical covariates should be also in covariates")
      }
      checkmate::assert_true(!("ss_auc" %in% covariates)) # special name
      checkmate::assert_true(!("dose" %in% covariates)) # special name
      checkmate::assert_class(system, "MultistateSystem")
      checkmate::assert_integerish(n_grid, lower = 10, len = 1)
      checkmate::assert_logical(prior_only, len = 1)
      checkmate::assert_logical(pk_only, len = 1)
      if (!is.null(pk_model)) {
        checkmate::assert_class(pk_model, "PKModel")
      }
      private$hazard_covariates <- covariates
      private$categorical <- categorical
      self$pk_model <- pk_model
      self$system <- system
      checkmate::assert_number(t_max, lower = 0)
      checkmate::assert_integerish(num_knots, lower = 3, upper = 20)
      self$set_knots(t_max, default_event_distribution(t_max), num_knots)
      private$n_grid <- n_grid
      self$prior_only <- prior_only
      self$pk_only <- pk_only
    },

    #' @description Set knot locations based on event times
    #'
    #' @param t_max Max time
    #' @param t_event Occurred event times
    #' @param num_knots Total number of knots. Includes the boundary knots.
    #' Number of spline basis functions will be \code{num_knots + 1}.
    #' @description The knots define how the spline basis functions are
    #' set.
    set_knots = function(t_max, t_event, num_knots) {
      checkmate::assert_integerish(num_knots, len = 1, lower = 3)
      checkmate::assert_number(t_max, lower = 0)
      checkmate::assert_numeric(t_event, lower = 0, upper = t_max, min.len = 3)
      tk <- place_internal_knots(t_max, num_knots - 2, t_event)
      knots <- c(0, tk, t_max)
      self$system$set_knots(knots)
    },

    #' @description Get knots
    get_knots = function() {
      self$system$get_knots()
    },

    #' @description Get max time
    get_tmax = function() {
      self$system$get_tmax()
    },

    #' @description Get names of the states
    get_states = function() {
      self$system$tm()$states
    },

    #' @description Is there a PK submodel?
    has_pk = function() {
      !is.null(self$pk_model)
    },

    #' @description Print the object
    #'
    #' @return nothing
    print = function() {
      covs <- self$covs()
      x1 <- paste0("A MultistateModel with:")
      x2 <- paste0(" - Hazard covariates: {", paste0(covs, collapse = ", "), "}")
      msg <- paste(x1, x2, "\n", sep = "\n")
      cat(msg)
      print(self$system)
      if (self$has_pk()) {
        print(self$pk_model)
      }
    },

    #' @description Get the hazard covariates (including steady-state exposure
    #' if PK model is included)
    covs = function() {
      x <- private$hazard_covariates
      if (self$has_pk()) {
        x <- c(x, "ss_auc")
      }
      unique(x)
    },

    #' @description Get all covariates that need to be given as data
    #' @param which Which subset to get?
    data_covs = function(which = NULL) {
      if (is.null(which)) {
        x <- private$hazard_covariates
        if (self$has_pk()) {
          x <- c(x, self$pk_model$covs())
        }
      } else {
        if (which == "haz") {
          x <- private$hazard_covariates
        } else {
          if (!self$has_pk()) {
            return(NULL)
          }
          if (which == "ka") {
            x <- self$pk_model$ka_covs()
          } else if (which == "CL") {
            x <- self$pk_model$CL_covs()
          } else if (which == "V2") {
            x <- self$pk_model$V2_covs()
          } else {
            if (!is.null(which)) {
              stop("invalid 'which' argument")
            }
          }
        }
      }

      unique(x)
    },

    #' @description Get names of categorical covariates
    categ_covs = function() {
      unique(private$categorical)
    },

    #' @description Simulate subject data, all covariates independently.
    #'
    #' @param N_subject Number of subjects.
    #' @param doses Possible doses.
    simulate_subjects = function(N_subject = 100, doses = c(15, 30, 60)) {
      checkmate::assert_numeric(doses, min.len = 1, lower = 0)

      # Generate covariates
      covs <- self$data_covs()
      categ <- self$categ_covs()
      idx_cat <- which(covs %in% categ)
      N_covs <- length(covs)
      A <- 60 + 10 * matrix(rnorm(N_subject * N_covs), N_subject, N_covs)

      # Discretize covariates
      for (idx in idx_cat) {
        A[, idx] <- as.numeric(A[, idx] > mean(A[, idx]))
      }

      # Create data frame
      df <- data.frame(A)
      colnames(df) <- covs
      df$subject_id <- sim_subject_ids(N_subject)
      n_groups <- length(doses)
      doses_vec <- doses[sample.int(n_groups, N_subject, replace = TRUE)]
      if (self$has_pk()) {
        df$dose <- doses_vec
      }
      if ("dose_amt" %in% covs) {
        df$dose_amt <- doses_vec
      }
      as_tibble(df)
    },

    #' @description Simulate data using the multistate model
    #'
    #' @param N_subject Number of subjects.
    #' @param beta_haz Covariate effects on each transition type.
    #' A matrix of shape \code{num_target_states} x \code{num_covs}.
    #' If \code{NULL}, a data frame of zeros is used.
    #' @param beta_pk Covariate effects on PK parameters. A named list with
    #' three elements, each being a vector. If any element is \code{NULL},
    #' a vector of zeros is used.
    #' @param w0 Baseline hazard rate for all transitions.
    #' @param w Spline weights. Matrix of shape \code{num_trans} x
    #' \code{num_weights}. If \code{NULL}, a matrix of zeros is used.
    #' @param num_doses Average number of doses taken by each subject. Only
    #' has effect if model as a PK submodel.
    #' @param subjects_df Subject data frame. If \code{NULL}, simulated using
    #' the \code{simulate_subjects} method.
    #' @return A \code{\link{JointData}} object.
    simulate_data = function(N_subject = 100, beta_haz = NULL,
                             beta_pk = NULL, w0 = 1e-3, w = NULL, num_doses = 10,
                             subjects_df = NULL) {
      H <- self$system$num_trans()
      if (is.null(subjects_df)) {
        subjects_df <- self$simulate_subjects(N_subject)
      } else {
        checkmate::assert_data_frame(subjects_df, nrows = N_subject)
      }

      checkmate::assert_numeric(w0, lower = 0)
      if (length(w0) > 1) {
        checkmate::assert_numeric(w0, len = H)
      } else {
        w0 <- rep(w0, H)
      }
      log_w0 <- log(w0)
      if (is.null(beta_haz)) {
        L <- length(self$target_states())
        K <- length(self$covs())
        beta_haz <- matrix(0, L, K)
      }
      sub_df_pk <- add_dosing_sim_opts(subjects_df, num_doses)
      pksim <- private$simulate_pk_data(sub_df_pk, beta_pk)
      pk_dat <- pksim$pk
      if (self$has_pk()) {
        subjects_df <- subjects_df |> dplyr::left_join(pk_dat, by = "subject_id")
      }
      path_df <- private$simulate_events(subjects_df, beta_haz, log_w0, w)
      link_df <- data.frame(
        path_id = seq_len(N_subject),
        subject_id = subjects_df$subject_id
      )
      link_df$rep_idx <- rep(1, N_subject)
      link_df$draw_idx <- rep(1, N_subject)
      pd <- PathData$new(
        subjects_df, path_df, link_df, self$system$tm(),
        colnames(subjects_df)
      )
      JointData$new(pd, pksim$dosing)
    },

    #' @description Get indices of states that are not source states
    #'
    #' @return integer
    target_states = function() {
      df <- self$system$tm()$states_df() |> dplyr::filter(!.data$source)
      df$state_idx
    }
  )
)

add_dosing_sim_opts <- function(df, num_doses) {
  checkmate::assert_integerish(num_doses, len = 1, lower = 4)
  df$num_doses <- num_doses - 3 + sample(5, size = nrow(df), replace = TRUE)
  df$num_ss_doses <- df$num_doses - 2
  df
}


# Create internal knots based on event time quantiles
place_internal_knots <- function(t_max, num_knots, t_event) {
  h <- 1 / (num_knots + 1)
  knots <- stats::quantile(t_event, probs = seq(0, 1, h))
  knots[2:(length(knots) - 1)]
}

# Normalize columns of matrix A
normalize_columns <- function(A) {
  K <- ncol(A)
  cn <- colnames(A)
  for (j in seq_len(K)) {
    A[, j] <- normalize(A[, j], cn[j])
  }
  A
}

# Normalize
normalize <- function(a, name) {
  sdd <- stats::sd(a)
  if (sdd == 0) {
    stop("error in normalization of ", name, ", zero variance")
  }
  x_norm <- (a - mean(a)) / sdd
  check_normalized_covariate(x_norm, name)
  x_norm
}
