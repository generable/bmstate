# Drop only first dimension
select_drop_first_dim_arr3 <- function(A, j) {
  A <- A[j, , , drop = FALSE]
  dim(A) <- dim(A)[-1]
  A
}

# Drop only first dimension
select_drop_first_dim_arr4 <- function(A, j) {
  A <- A[j, , , , drop = FALSE]
  dim(A) <- dim(A)[-1]
  A
}

# Utility
rv <- function(fit, name) {
  posterior::as_draws_rvars(fit$draws(name))[[name]]
}
