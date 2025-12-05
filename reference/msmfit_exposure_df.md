# Compute exposure and return as rvar in df

Compute exposure and return as rvar in df

## Usage

``` r
msmfit_exposure_df(fit, oos = FALSE, data = NULL)
```

## Arguments

- fit:

  A
  [`MultistateModelFit`](https://generable.github.io/bmstate/reference/MultistateModelFit.md)
  object

- oos:

  Out-of-sample mode? If `FALSE`, the possible subject-specific fitted
  parameters are used. If `TRUE`, acting as if the subjects are new.

- data:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  object. If `NULL`, the data used to fit the model is used.
