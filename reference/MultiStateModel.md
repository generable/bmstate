# Main model class

Class that represents a multistate model.

## Public fields

- `system`:

  A
  [`MultistateSystem`](https://generable.github.io/bmstate/reference/MultistateSystem.md)

- `pk_model`:

  A
  [`PKModel`](https://generable.github.io/bmstate/reference/PKModel.md)
  or NULL.

- `prior_only`:

  Should the model ignore likelihood?

- `pk_only`:

  Should the model ignore the entire hazard model part?

## Methods

### Public methods

- [`MultistateModel$get_normalizers()`](#method-MultistateModel-get_normalizers)

- [`MultistateModel$get_n_grid()`](#method-MultistateModel-get_n_grid)

- [`MultistateModel$set_normalizers()`](#method-MultistateModel-set_normalizers)

- [`MultistateModel$get_auc_normalizers()`](#method-MultistateModel-get_auc_normalizers)

- [`MultistateModel$set_auc_normalizers()`](#method-MultistateModel-set_auc_normalizers)

- [`MultistateModel$get_prior_mean_h0()`](#method-MultistateModel-get_prior_mean_h0)

- [`MultistateModel$set_prior_mean_h0()`](#method-MultistateModel-set_prior_mean_h0)

- [`MultistateModel$set_prior_mean_h0_data()`](#method-MultistateModel-set_prior_mean_h0_data)

- [`MultistateModel$new()`](#method-MultistateModel-new)

- [`MultistateModel$set_knots()`](#method-MultistateModel-set_knots)

- [`MultistateModel$get_knots()`](#method-MultistateModel-get_knots)

- [`MultistateModel$get_tmax()`](#method-MultistateModel-get_tmax)

- [`MultistateModel$get_states()`](#method-MultistateModel-get_states)

- [`MultistateModel$has_pk()`](#method-MultistateModel-has_pk)

- [`MultistateModel$print()`](#method-MultistateModel-print)

- [`MultistateModel$covs()`](#method-MultistateModel-covs)

- [`MultistateModel$data_covs()`](#method-MultistateModel-data_covs)

- [`MultistateModel$categ_covs()`](#method-MultistateModel-categ_covs)

- [`MultistateModel$simulate_subjects()`](#method-MultistateModel-simulate_subjects)

- [`MultistateModel$simulate_data()`](#method-MultistateModel-simulate_data)

- [`MultistateModel$target_states()`](#method-MultistateModel-target_states)

- [`MultistateModel$clone()`](#method-MultistateModel-clone)

------------------------------------------------------------------------

### Method `get_normalizers()`

Get normalization constants for each variable

#### Usage

    MultistateModel$get_normalizers()

#### Returns

list

------------------------------------------------------------------------

### Method `get_n_grid()`

Get number of grid points used for integration.

#### Usage

    MultistateModel$get_n_grid()

#### Returns

An integer

------------------------------------------------------------------------

### Method `set_normalizers()`

Set normalization constant for each variable (side effect)

#### Usage

    MultistateModel$set_normalizers(data)

#### Arguments

- `data`:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  object

------------------------------------------------------------------------

### Method `get_auc_normalizers()`

Get normalization constants for AUC (PK)

#### Usage

    MultistateModel$get_auc_normalizers()

#### Returns

list

------------------------------------------------------------------------

### Method `set_auc_normalizers()`

Set normalization constants for AUC (side effect)

#### Usage

    MultistateModel$set_auc_normalizers(loc = 0, scale = 1)

#### Arguments

- `loc`:

  Location

- `scale`:

  Scale

------------------------------------------------------------------------

### Method `get_prior_mean_h0()`

Get assumed prior mean baseline hazard rates.

#### Usage

    MultistateModel$get_prior_mean_h0()

#### Returns

Numeric vector with length equal to number of transitions

------------------------------------------------------------------------

### Method `set_prior_mean_h0()`

Set assumed prior mean baseline hazard rates (side effect).

#### Usage

    MultistateModel$set_prior_mean_h0(mean_h0)

#### Arguments

- `mean_h0`:

  Numeric vector with length equal to number of transitions

------------------------------------------------------------------------

### Method `set_prior_mean_h0_data()`

Set assumed prior mean baseline hazard rates (side effect) based on
average hazards in data.

#### Usage

    MultistateModel$set_prior_mean_h0_data(data)

#### Arguments

- `data`:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  or
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object.

------------------------------------------------------------------------

### Method `new()`

Create model

#### Usage

    MultistateModel$new(
      system,
      covariates = NULL,
      pk_model = NULL,
      t_max = 1000,
      num_knots = 5,
      categorical = NULL,
      n_grid = 1000,
      prior_only = FALSE,
      pk_only = FALSE
    )

#### Arguments

- `system`:

  A
  [`MultistateSystem`](https://generable.github.io/bmstate/reference/MultistateSystem.md)

- `covariates`:

  The names of the hazard covariates (excluding possible exposure
  estimated from PK model). Do not use reserved names `ss_auc` or
  `dose`.

- `pk_model`:

  A
  [`PKModel`](https://generable.github.io/bmstate/reference/PKModel.md)
  or NULL.

- `t_max`:

  Max time.

- `num_knots`:

  Total number of spline knots.

- `categorical`:

  Names of covariates that are binary. This only has an effect when
  simulating data. When fitting a model, all covariates are treated as
  continuous, so you should use a binary encoding for categories if
  there is more than two.

- `n_grid`:

  Number of time discretization points for integrating

- `prior_only`:

  Should the model ignore likelihood?

- `pk_only`:

  Should the model ignore the entire hazard model part? hazards.

------------------------------------------------------------------------

### Method `set_knots()`

Set knot locations based on event times

The knots define how the spline basis functions are set.

#### Usage

    MultistateModel$set_knots(t_max, t_event, num_knots)

#### Arguments

- `t_max`:

  Max time

- `t_event`:

  Occurred event times

- `num_knots`:

  Total number of knots. Includes the boundary knots. Number of spline
  basis functions will be `num_knots + 1`.

------------------------------------------------------------------------

### Method `get_knots()`

Get knots

#### Usage

    MultistateModel$get_knots()

------------------------------------------------------------------------

### Method `get_tmax()`

Get max time

#### Usage

    MultistateModel$get_tmax()

------------------------------------------------------------------------

### Method `get_states()`

Get names of the states

#### Usage

    MultistateModel$get_states()

------------------------------------------------------------------------

### Method `has_pk()`

Is there a PK submodel?

#### Usage

    MultistateModel$has_pk()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the object

#### Usage

    MultistateModel$print()

#### Returns

nothing

------------------------------------------------------------------------

### Method `covs()`

Get the hazard covariates (including steady-state exposure if PK model
is included)

#### Usage

    MultistateModel$covs()

------------------------------------------------------------------------

### Method `data_covs()`

Get all covariates that need to be given as data

#### Usage

    MultistateModel$data_covs(which = NULL)

#### Arguments

- `which`:

  Which subset to get?

------------------------------------------------------------------------

### Method `categ_covs()`

Get names of categorical covariates

#### Usage

    MultistateModel$categ_covs()

------------------------------------------------------------------------

### Method `simulate_subjects()`

Simulate subject data, all covariates independently.

#### Usage

    MultistateModel$simulate_subjects(N_subject = 100, doses = c(15, 30, 60))

#### Arguments

- `N_subject`:

  Number of subjects.

- `doses`:

  Possible doses.

------------------------------------------------------------------------

### Method `simulate_data()`

Simulate data using the multistate model

#### Usage

    MultistateModel$simulate_data(
      N_subject = 100,
      beta_haz = NULL,
      beta_pk = NULL,
      w0 = 0.001,
      w = NULL,
      num_doses = 10,
      subjects_df = NULL
    )

#### Arguments

- `N_subject`:

  Number of subjects.

- `beta_haz`:

  Covariate effects on each transition type. A matrix of shape
  `num_target_states` x `num_covs`. If `NULL`, a data frame of zeros is
  used.

- `beta_pk`:

  Covariate effects on PK parameters. A named list with three elements,
  each being a vector. If any element is `NULL`, a vector of zeros is
  used.

- `w0`:

  Baseline hazard rate for all transitions.

- `w`:

  Spline weights. Matrix of shape `num_trans` x `num_weights`. If
  `NULL`, a matrix of zeros is used.

- `num_doses`:

  Average number of doses taken by each subject. Only has effect if
  model as a PK submodel.

- `subjects_df`:

  Subject data frame. If `NULL`, simulated using the `simulate_subjects`
  method.

#### Returns

A
[`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
object.

------------------------------------------------------------------------

### Method `target_states()`

Get indices of states that are not source states

#### Usage

    MultistateModel$target_states()

#### Returns

integer

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultistateModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
