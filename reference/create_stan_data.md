# Creating Stan data list

Creating Stan data list

## Usage

``` r
create_stan_data(model, data)
```

## Arguments

- model:

  A
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)
  object.

- data:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  or
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object.

## Value

A list of data for Stan.

## See also

Other Stan-related functions:
[`create_stan_model()`](https://generable.github.io/bmstate/reference/create_stan_model.md),
[`default_stan_filepath()`](https://generable.github.io/bmstate/reference/default_stan_filepath.md),
[`ensure_exposed_stan_functions()`](https://generable.github.io/bmstate/reference/ensure_exposed_stan_functions.md),
[`fit_stan()`](https://generable.github.io/bmstate/reference/fit_stan.md),
[`plot_stan_data_integral()`](https://generable.github.io/bmstate/reference/plot_stan_data_integral.md),
[`plot_stan_data_matrix()`](https://generable.github.io/bmstate/reference/plot_stan_data_matrix.md)
