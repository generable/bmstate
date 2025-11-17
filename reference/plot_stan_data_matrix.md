# Visualize parts of the 'Stan' data

Visualize parts of the 'Stan' data

## Usage

``` r
plot_stan_data_matrix(model, stan_data, name, subject_idx)
```

## Arguments

- model:

  A
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md).

- stan_data:

  A list of data for Stan. Created using
  [`create_stan_data`](https://generable.github.io/bmstate/reference/create_stan_data.md).

- name:

  Name of Stan data field to plot. Currently must be one of
  `"transition", "at_risk"`.

- subject_idx:

  Index of subject whose data rows to plot.

## See also

Other Stan-related functions:
[`create_stan_data()`](https://generable.github.io/bmstate/reference/create_stan_data.md),
[`create_stan_model()`](https://generable.github.io/bmstate/reference/create_stan_model.md),
[`default_stan_filepath()`](https://generable.github.io/bmstate/reference/default_stan_filepath.md),
[`ensure_exposed_stan_functions()`](https://generable.github.io/bmstate/reference/ensure_exposed_stan_functions.md),
[`fit_stan()`](https://generable.github.io/bmstate/reference/fit_stan.md),
[`plot_stan_data_integral()`](https://generable.github.io/bmstate/reference/plot_stan_data_integral.md)
