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

# Default event distribution
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
    dplyr::filter(.data$state %in% term_state_inds, is_event == 1) |>
    dplyr::group_by(.data$path_id) |>
    summarise(term_time = min(.data$time, na.rm = T)) |>
    dplyr::ungroup()
  no_terms <- df |>
    dplyr::anti_join(term_events, by = "path_id")
  with_terms <- df |>
    inner_join(term_events, by = c("path_id")) |>
    dplyr::filter(.data$time <= term_time) |>
    dplyr::select(-"term_time")
  no_terms |>
    dplyr::bind_rows(with_terms)
}
