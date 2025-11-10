# Solve state occupancy probabilities for each subject and draw in 'MultistateModelFit'

Solve state occupancy probabilities for each subject and draw in
'MultistateModelFit'

## Usage

``` r
p_state_occupancy(
  fit,
  oos = FALSE,
  t_start = 0,
  t_out = NULL,
  data = NULL,
  ...
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

  Start time. The results are transition probabilities to each state,
  evaluated at each time point of `t_out`, given that the that the
  subject is in the state they are in `data` at `t_start`.

- t_out:

  End times (vector)

- data:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  object. If `NULL`, the data used to fit the model is used.

- ...:

  Arguments passed to
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Value

A data frame with columns

- `subject_id` (character)

- `time` (numeric)

- `state` (character)

- `prob` state occupancy probability for given subject at given time

Note that `prob` is an `rvar` over all parameter draws.
