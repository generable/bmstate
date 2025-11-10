# The 'bmstate' package.

Bayesian multistate modeling of right-censored time-to-event data.

## Transition matrix

For any analysis with this package, you will need a
[`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)
which specifies the states and possible transitions between them. This
can be created directly with `TransitionMatrix$new()` or using built-in
functions for common multistate systems such as
[`transmat_illnessdeath`](https://generable.github.io/bmstate/reference/transmat.md).

## Data format

See the function
[`df_to_pathdata`](https://generable.github.io/bmstate/reference/df_to_pathdata.md)
for how to create a
[`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
object from your data frame and a given
[`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md).

## Creating a model

You can create a model using
[`create_msm`](https://generable.github.io/bmstate/reference/create_msm.md).
A model is represented by a
[`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)
object, which has a
[`MultistateSystem`](https://generable.github.io/bmstate/reference/MultistateSystem.md),
which has a
[`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md).

## Fitting a model

You can fit a model with 'Stan' using
[`fit_stan`](https://generable.github.io/bmstate/reference/fit_stan.md).

## Prediction with a fitted model

- For predicting state visit risks over a time period using state path
  simulation, see the functions
  [`generate_paths`](https://generable.github.io/bmstate/reference/generate_paths.md),
  [`p_state_visit`](https://generable.github.io/bmstate/reference/p_state_visit.md),
  and
  [`p_state_visit_per_subject`](https://generable.github.io/bmstate/reference/p_state_visit_per_subject.md).

- For analytically solving state occupancy probabilities after some
  time, see
  [`p_state_occupancy`](https://generable.github.io/bmstate/reference/p_state_occupancy.md).

## Vignettes

See the vignettes for more information.

## See also

Useful links:

- <https://generable.github.io/bmstate/>

## Author

Juho Timonen and Jacqueline Buros
