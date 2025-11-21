# Look for potential covariates that affect transitions

Look for potential covariates that affect transitions

## Usage

``` r
potential_covariates(pd, possible = NULL, ...)
```

## Arguments

- pd:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object.

- possible:

  Possible covariates to look for (character vector). If `NULL`
  (default) all covariates of the object are tested.

- ...:

  Arguments passed to `fit_coxph()`.

## Value

A `data.frame` with columns `pval`, `target_state`, `covariate`.
