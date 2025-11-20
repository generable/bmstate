# Fit a model using 'Stan'

*NOTE:* This function has a side effect of setting covariate
normalizers, prior assumed mean baseline hazard, and concentration upper
bound (PK) based on data.

## Usage

``` r
fit_stan(
  model,
  data,
  set_auc_normalizers = TRUE,
  filepath = NULL,
  return_stanfit = FALSE,
  max_conc_factor = 100,
  method = "sample",
  ...
)
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

- set_auc_normalizers:

  Set AUC normalization based on SS doses.

- filepath:

  Deprecated.

- return_stanfit:

  Return also the raw 'Stan' fit object?

- max_conc_factor:

  Factor to multiply observed max concentration by to get concentration
  upper bound.

- method:

  Must be one of `"sample"` (default), `"pathfinder"` or `"optimize"`.

- ...:

  Arguments passed to the `sample`, `pathfinder` or `optimize` method of
  the 'CmdStanR' model.

## Value

A
[`MultistateModelFit`](https://generable.github.io/bmstate/reference/MultistateModelFit.md)
object.

## See also

Other Stan-related functions:
[`create_stan_data()`](https://generable.github.io/bmstate/reference/create_stan_data.md),
[`create_stan_model()`](https://generable.github.io/bmstate/reference/create_stan_model.md),
[`default_stan_filepath()`](https://generable.github.io/bmstate/reference/default_stan_filepath.md),
[`ensure_exposed_stan_functions()`](https://generable.github.io/bmstate/reference/ensure_exposed_stan_functions.md),
[`plot_stan_data_integral()`](https://generable.github.io/bmstate/reference/plot_stan_data_integral.md),
[`plot_stan_data_matrix()`](https://generable.github.io/bmstate/reference/plot_stan_data_matrix.md)
