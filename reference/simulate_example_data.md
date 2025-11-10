# Simulate example data

Simulate example data

## Usage

``` r
simulate_example_data(
  N = 10,
  beta_haz = NULL,
  beta_pk = NULL,
  w0 = 0.001,
  tm = NULL,
  ...
)
```

## Arguments

- N:

  number of subjects

- beta_haz:

  see the documentation of
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)

- beta_pk:

  see the documentation of
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)

- w0:

  see the documentation of
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)

- tm:

  A
  [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)

- ...:

  arguments passed to
  [`create_msm`](https://generable.github.io/bmstate/reference/create_msm.md)
