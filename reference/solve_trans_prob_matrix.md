# Solve the transition probability matrix for each given output time

Solve the transition probability matrix for each given output time

## Usage

``` r
solve_trans_prob_matrix(
  system,
  t_out,
  log_w0,
  w = NULL,
  log_m = NULL,
  t_start = 0,
  ...
)
```

## Arguments

- system:

  A
  [`MultistateSystem`](https://generable.github.io/bmstate/reference/MultistateSystem.md)

- t_out:

  End times (vector)

- log_w0:

  A vector of length `n_trans`

- w:

  An array of shape `n_trans` x `n_weights`

- log_m:

  A vector of length `n_trans`

- t_start:

  Initial time

- ...:

  Arguments passed to
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Value

An array `P` where `P[k,i,j]` is the probability that the system will be
in state `j` at the `k`th time point of `t_out` given that it is in
state `i` at time `t_start`.
