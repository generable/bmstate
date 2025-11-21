# Create a multistate model

Create a multistate model

## Usage

``` r
create_msm(tm, hazard_covs = NULL, pk_covs = NULL, categ_covs = NULL, ...)
```

## Arguments

- tm:

  A
  [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md).
  See
  [`transmat`](https://generable.github.io/bmstate/reference/transmat.md)
  for how to create common transition matrices.

- hazard_covs:

  Covariates that affect the hazard. A character vector. The name
  `"dose_amt"` is special if simulating data using the model.

- pk_covs:

  Covariates that affect the PK parameters. A list with elements `ka`
  `CL`, and `V2`. If `NULL`, a PK model will not be created.

- categ_covs:

  Names of covariates that are binary. This only has an effect when
  simulating data. When fitting a model, all covariates are treated as
  continuous, so you should use a binary encoding for categories if
  there is more than two.

- ...:

  Arguments passed to
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)
  init

## Value

A
[`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)
object.
