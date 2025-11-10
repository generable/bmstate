# Extract and reshape draws of instant hazard related parameters

Extract and reshape draws of instant hazard related parameters

## Usage

``` r
msmfit_inst_hazard_param_draws(fit, oos = FALSE, data = NULL)
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

## Value

a list with elements `log_m`, `log_w0`, `w`, each of which is an array
where the first dimension is the number of subjects times the number of
draws
