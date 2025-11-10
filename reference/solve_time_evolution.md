# Solve the Kolmogorov forward equation of a Markovian multistate system

Solve the Kolmogorov forward equation of a Markovian multistate system

## Usage

``` r
solve_time_evolution(system, t, log_w0, w = NULL, log_m = NULL, ...)
```

## Arguments

- system:

  A
  [`MultistateSystem`](https://generable.github.io/bmstate/reference/MultistateSystem.md)

- t:

  A time grid

- log_w0:

  A vector of length `n_trans`

- w:

  An array of shape `n_trans` x `n_weights`

- log_m:

  A vector of length `n_trans`

- ...:

  Arguments passed to
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Value

Value given by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).
