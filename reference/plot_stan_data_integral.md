# Visualize numerical integration of hazard

Visualize numerical integration of hazard

## Usage

``` r
plot_stan_data_integral(model, stan_data, iidx, h0 = 0.001, w = NULL)
```

## Arguments

- model:

  A
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md).

- stan_data:

  A list of data for Stan. Created using
  [`create_stan_data`](https://generable.github.io/bmstate/reference/create_stan_data.md).

- iidx:

  Interval index

- h0:

  Hazard function average value

- w:

  Hazard function spline weights. Randomly generated if not given.

## See also

Other Stan-related functions:
[`create_stan_model()`](https://generable.github.io/bmstate/reference/create_stan_model.md),
[`default_stan_filepath()`](https://generable.github.io/bmstate/reference/default_stan_filepath.md),
[`ensure_exposed_stan_functions()`](https://generable.github.io/bmstate/reference/ensure_exposed_stan_functions.md),
[`fit_stan()`](https://generable.github.io/bmstate/reference/fit_stan.md),
[`plot_stan_data_matrix()`](https://generable.github.io/bmstate/reference/plot_stan_data_matrix.md)
