# Create the main 'Stan' model

Uses the Stan code at at
[`default_stan_filepath`](https://generable.github.io/bmstate/reference/default_stan_filepath.md)
or a custom file if it has been set using
`options(bmstate_stan_file = ...)`.

## Usage

``` r
create_stan_model(filepath = NULL, ...)
```

## Arguments

- filepath:

  Deprecated.

- ...:

  Arguments passed to 'CmdStanR'

## See also

Other Stan-related functions:
[`create_stan_data()`](https://generable.github.io/bmstate/reference/create_stan_data.md),
[`default_stan_filepath()`](https://generable.github.io/bmstate/reference/default_stan_filepath.md),
[`ensure_exposed_stan_functions()`](https://generable.github.io/bmstate/reference/ensure_exposed_stan_functions.md),
[`fit_stan()`](https://generable.github.io/bmstate/reference/fit_stan.md),
[`plot_stan_data_integral()`](https://generable.github.io/bmstate/reference/plot_stan_data_integral.md),
[`plot_stan_data_matrix()`](https://generable.github.io/bmstate/reference/plot_stan_data_matrix.md)
