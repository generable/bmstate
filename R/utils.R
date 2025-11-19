# Check for unusual normalized covariate values (will cause issues for hazard)
check_normalized_covariate <- function(x_norm) {
  mabs <- max(abs(x_norm))
  msg <- paste0(
    "a normalized covariate has maximum absolute value ",
    round(mabs, 5), ", are you sure normalization of covariates is correct?"
  )
  if (mabs > 10) {
    warning(msg)
    return(FALSE)
  }
  TRUE
}


# Access one dimension and drop only that dimension
access_one_dim <- function(x, dim_index, value) {
  dims <- rep(list(quote(expr = )), length(dim(x)))
  dims[[dim_index]] <- value
  result <- eval(as.call(c(quote(`[`), quote(x), dims, drop = FALSE)))

  # Remove the singleton dimension manually
  dim(result) <- dim(result)[-dim_index]
  result
}

# Utility
rv <- function(fit, name) {
  posterior::as_draws_rvars(fit$draws(name))[[name]]
}

# Helper
create_rv_list <- function(stan_fit, names) {
  out <- list()
  j <- 0
  names_out <- NULL
  for (name in names) {
    tryCatch(
      {
        rv <- rv(stan_fit, name)
        j <- j + 1
        out[[j]] <- rv
        names_out <- c(names_out, name)
      },
      error = function(e) {
      }
    )
  }
  names(out) <- names_out
  out
}

# Rvar to an rvar with a single draw corresponding to the mean of original draws
rvar_to_mean_rvar <- function(rv) {
  mrv <- mean(rv)
  D <- dim(mrv)
  if (is.null(D)) {
    L <- length(mrv)
    out <- posterior::rvar(array(mrv, dim = c(1, L)), dim = L)
  } else {
    out <- posterior::rvar(array(mrv, dim = c(1, D)), dim = D)
  }
  out
}


# Ordered by value
find_row_and_col_of_positive_vals <- function(mat) {
  positions <- which(mat > 0, arr.ind = TRUE)
  values <- mat[positions]
  out <- positions[order(values), ]
  if (is.null(dim(out))) {
    out <- as.matrix(t(out))
  }
  out
}

# Default event time distribution
default_event_distribution <- function(t_max) {
  x <- seq(0, t_max, length.out = 100)
  max(x) * x^2 / (max(x^2))
}

# Subject ids for simulated data
sim_subject_ids <- function(N) {
  paste0("sim-", formatC(seq_len(N), 5, flag = "0"))
}

# Truncate path data frame
truncate_after_terminal_events <- function(df, term_state_inds) {
  term_events <- df |>
    dplyr::filter(.data$state %in% term_state_inds, .data$trans_idx > 0) |>
    dplyr::group_by(.data$path_id) |>
    summarise(term_time = min(.data$time, na.rm = T)) |>
    dplyr::ungroup()
  no_terms <- df |>
    dplyr::anti_join(term_events, by = "path_id")
  with_terms <- df |>
    inner_join(term_events, by = c("path_id")) |>
    dplyr::filter(.data$time <= .data$term_time) |>
    dplyr::select(-"term_time")
  no_terms |>
    dplyr::bind_rows(with_terms)
}

# Checks to be made that make sure that the model settings
# make sense for given data
prefit_checks <- function(model, data) {
  pd <- data$paths
  lens <- pd$as_transitions() |>
    dplyr::mutate(time_len = .data$time - .data$time_prev) |>
    dplyr::pull(.data$time_len)

  # Numbers that can lead to message, warning, or error if bad
  tmax <- model$get_tmax()
  min_len <- min(lens)
  delta_grid <- tmax / model$get_n_grid()
  max_time <- max(pd$get_path_df()$time)

  # Define possible messages
  msg <- paste0(
    "Shortest time interval (", min_len,
    ") is smaller than delta_grid (", delta_grid,
    "). Consider increasing n_grid or decreasing t_max of the model."
  )
  msg_warning <- paste0(
    "Set model t_max is more than twice as large as ",
    "largest observed time in the data. Are you sure you want to do this? ",
    "If not, set smaller t_max."
  )
  msg_error <- paste0(
    "Model t_max (", tmax,
    ") is smaller than max observed time in data (",
    max_time, "). Increase t_max."
  )

  # Check
  if (min_len < delta_grid) {
    message(msg)
  }
  if (tmax > 2 * max_time) {
    warning(msg_warning)
  }
  if (tmax < max_time) {
    stop(msg_error)
  }
  TRUE
}
