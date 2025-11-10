# Path generation for 'MultistateModelFit'

Generates paths for each subject in the data, starting from the state
that the subject is at the given start time.

## Usage

``` r
generate_paths(
  fit,
  oos = FALSE,
  t_start = 0,
  t_max = NULL,
  n_rep = 10,
  data = NULL
)
```

## Arguments

- fit:

  A
  [`MultistateModelFit`](https://generable.github.io/bmstate/reference/MultistateModelFit.md)
  object

- oos:

  Out-of-sample mode? If `FALSE`, the possible subject-specific fitted
  parameters are used. If `TRUE`, acting as if the subjects are new.

- t_start:

  Start time. A single value or a vector with length equal to number of
  subjects.

- t_max:

  Max time. If `NULL`, the max time of the model is used.

- n_rep:

  Number of repeats per draw.

- data:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  object. If `NULL`, the data used to fit the model is used.

## Value

A
[`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
object.
