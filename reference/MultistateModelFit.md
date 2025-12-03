# Minimal fit class

Class that represents a multistate model with certain parameters or
parameter draws.

## Public fields

- `model`:

  The
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)

- `data`:

  A
  [`JointData`](https://generable.github.io/bmstate/reference/JointData.md)
  object

- `info`:

  Fit info.

## Methods

### Public methods

- [`MultistateModelFit$new()`](#method-MultistateModelFit-new)

- [`MultistateModelFit$is_pk_only()`](#method-MultistateModelFit-is_pk_only)

- [`MultistateModelFit$assert_hazard_fit()`](#method-MultistateModelFit-assert_hazard_fit)

- [`MultistateModelFit$mean_fit()`](#method-MultistateModelFit-mean_fit)

- [`MultistateModelFit$is_point_estimate()`](#method-MultistateModelFit-is_point_estimate)

- [`MultistateModelFit$get_data()`](#method-MultistateModelFit-get_data)

- [`MultistateModelFit$draws_names()`](#method-MultistateModelFit-draws_names)

- [`MultistateModelFit$get_draws()`](#method-MultistateModelFit-get_draws)

- [`MultistateModelFit$get_draws_of()`](#method-MultistateModelFit-get_draws_of)

- [`MultistateModelFit$print()`](#method-MultistateModelFit-print)

- [`MultistateModelFit$plot_basisfun()`](#method-MultistateModelFit-plot_basisfun)

- [`MultistateModelFit$simulate_pk()`](#method-MultistateModelFit-simulate_pk)

- [`MultistateModelFit$plot_pk()`](#method-MultistateModelFit-plot_pk)

- [`MultistateModelFit$plot_h0()`](#method-MultistateModelFit-plot_h0)

- [`MultistateModelFit$h0_dist()`](#method-MultistateModelFit-h0_dist)

- [`MultistateModelFit$covariate_effects()`](#method-MultistateModelFit-covariate_effects)

- [`MultistateModelFit$log_z_pars()`](#method-MultistateModelFit-log_z_pars)

- [`MultistateModelFit$num_draws()`](#method-MultistateModelFit-num_draws)

- [`MultistateModelFit$clone()`](#method-MultistateModelFit-clone)

------------------------------------------------------------------------

### Method `new()`

Create model fit object

#### Usage

    MultistateModelFit$new(data, stan_data, model, draws, info = NULL)

#### Arguments

- `data`:

  Data used to create the fit.

- `stan_data`:

  The used 'Stan' data list.

- `model`:

  A
  [`MultistateModel`](https://generable.github.io/bmstate/reference/MultiStateModel.md)

- `draws`:

  A named list of rvars.

- `info`:

  Fit info.

------------------------------------------------------------------------

### Method `is_pk_only()`

Is this a PK-only fit?

#### Usage

    MultistateModelFit$is_pk_only()

#### Returns

logical

------------------------------------------------------------------------

### Method `assert_hazard_fit()`

Require that fit is not PK-only

#### Usage

    MultistateModelFit$assert_hazard_fit()

------------------------------------------------------------------------

### Method `mean_fit()`

Create a reduced version with only a single draw, corresponding to the
mean of original draws.

#### Usage

    MultistateModelFit$mean_fit()

#### Returns

A new `MultistateModelFit` object.

------------------------------------------------------------------------

### Method `is_point_estimate()`

Check if fit is a point estimate

#### Usage

    MultistateModelFit$is_point_estimate()

#### Returns

a logical value

------------------------------------------------------------------------

### Method `get_data()`

Extract Stan data list

#### Usage

    MultistateModelFit$get_data()

------------------------------------------------------------------------

### Method `draws_names()`

Names of the draws list

#### Usage

    MultistateModelFit$draws_names()

------------------------------------------------------------------------

### Method `get_draws()`

Extract draws as `rvar`s

#### Usage

    MultistateModelFit$get_draws(name = NULL)

#### Arguments

- `name`:

  Param/quantity name

------------------------------------------------------------------------

### Method `get_draws_of()`

Draws in a raw array with same shape as Stan variable

#### Usage

    MultistateModelFit$get_draws_of(name)

#### Arguments

- `name`:

  Param/quantity name of `x`

#### Returns

Array with dimension `c(ndraws(x), dim(x))`

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the object

#### Usage

    MultistateModelFit$print()

#### Returns

nothing

------------------------------------------------------------------------

### Method `plot_basisfun()`

Plot used basis functions (grid)

#### Usage

    MultistateModelFit$plot_basisfun()

------------------------------------------------------------------------

### Method `simulate_pk()`

Simulate PK dynamics using fitted params.

#### Usage

    MultistateModelFit$simulate_pk(
      oos = FALSE,
      data = NULL,
      L = 100,
      timescale = 24,
      n_prev = 3
    )

#### Arguments

- `oos`:

  Out-of-sample subjects?

- `data`:

  Data for which to predict the concentration. If `NULL`, training data
  is used.

- `L`:

  number of grid points for each subject

- `timescale`:

  scale of time

- `n_prev`:

  number of previous doses to show fit for

------------------------------------------------------------------------

### Method `plot_pk()`

Plot PK fit.

#### Usage

    MultistateModelFit$plot_pk(
      max_num_subjects = 12,
      oos = FALSE,
      data = NULL,
      L = 100,
      timescale = 24,
      n_prev = 3,
      ci_alpha = 0.9,
      subject_ids = NULL
    )

#### Arguments

- `max_num_subjects`:

  Max number of subjects to show.

- `oos`:

  Out-of-sample subjects?

- `data`:

  Data for which to predict the concentration. If `NULL`, training data
  is used.

- `L`:

  number of grid points for each subject

- `timescale`:

  scale of time

- `n_prev`:

  number of previous doses to show fit for

- `ci_alpha`:

  Width of central credible interval.

- `subject_ids`:

  Which subjects to plot?

------------------------------------------------------------------------

### Method `plot_h0()`

Plot baseline hazard distribution

#### Usage

    MultistateModelFit$plot_h0(t = NULL, ci_alpha = 0.95)

#### Arguments

- `t`:

  times where to evaluate the baseline hazards

- `ci_alpha`:

  width of credible interval

------------------------------------------------------------------------

### Method `h0_dist()`

Baseline hazard distribution

#### Usage

    MultistateModelFit$h0_dist(t = NULL, ci_alpha = 0.95)

#### Arguments

- `t`:

  times where to evaluate the baseline hazards

- `ci_alpha`:

  width of credible interval

------------------------------------------------------------------------

### Method `covariate_effects()`

Extract covariate effects

Currently not implemented for models that have a PK submodel.

#### Usage

    MultistateModelFit$covariate_effects()

#### Returns

A data frame which has columns

- `covariate` Name of the covariate

- `beta` The covariate effect parameter estimate (`rvar`). NOTE: this is
  the regression coefficient for *normalized* covariates.

- `target_state` Name of the target state. The corresponding beta is the
  covariate effect on all transitions that end in this target state.

------------------------------------------------------------------------

### Method `log_z_pars()`

Full names of parameters that start with `log_z_`.

#### Usage

    MultistateModelFit$log_z_pars()

------------------------------------------------------------------------

### Method `num_draws()`

Get number of draws

#### Usage

    MultistateModelFit$num_draws()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultistateModelFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
