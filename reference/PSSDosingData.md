# Partially steady-state dosing data class (R6 class)

Contains information about taken doses and dose times.

## Super class

[`bmstate::DosingData`](https://generable.github.io/bmstate/reference/DosingData.md)
-\> `PSSDosingData`

## Public fields

- `doses`:

  Dose amounts (list with length equal to number of subjects).
  Corresponds to doses after the steady state.

- `times`:

  Dose times (list with length equal to number of subjects). Corresponds
  to doses after the steady state. First time here is the end of
  steady-state assumption time range.

- `dose_ss`:

  Steady-state dose.

- `tau_ss`:

  Steady-state dosing interval.

## Methods

### Public methods

- [`PSSDosingData$new()`](#method-PSSDosingData-new)

- [`PSSDosingData$as_data_frame()`](#method-PSSDosingData-as_data_frame)

- [`PSSDosingData$simulate_pk()`](#method-PSSDosingData-simulate_pk)

- [`PSSDosingData$plot()`](#method-PSSDosingData-plot)

- [`PSSDosingData$filter()`](#method-PSSDosingData-filter)

- [`PSSDosingData$clone()`](#method-PSSDosingData-clone)

Inherited methods

- [`bmstate::DosingData$num_subjects()`](https://generable.github.io/bmstate/reference/DosingData.html#method-num_subjects)
- [`bmstate::DosingData$print()`](https://generable.github.io/bmstate/reference/DosingData.html#method-print)
- [`bmstate::DosingData$set_sub_ids()`](https://generable.github.io/bmstate/reference/DosingData.html#method-set_sub_ids)

------------------------------------------------------------------------

### Method `new()`

Initialize

#### Usage

    PSSDosingData$new(subject_ids, doses, times, dose_ss = NULL, tau_ss = 24)

#### Arguments

- `subject_ids`:

  A character vector

- `doses`:

  Dose amounts (list with length equal to number of subjects).
  Corresponds to doses after the steady state.

- `times`:

  Dose times (list with length equal to number of subjects). Corresponds
  to doses after the steady state. First time here is the end of
  steady-state assumption time range.

- `dose_ss`:

  Steady-state dose.

- `tau_ss`:

  Steady-state dosing interval.

------------------------------------------------------------------------

### Method `as_data_frame()`

As data frame

#### Usage

    PSSDosingData$as_data_frame()

------------------------------------------------------------------------

### Method `simulate_pk()`

Simulate PK dynamics

#### Usage

    PSSDosingData$simulate_pk(t, theta, MAX_CONC)

#### Arguments

- `t`:

  A vector of output times for each subject (a list).

- `theta`:

  A matrix of parameters.

- `MAX_CONC`:

  concentration upper bound

#### Returns

a `data.frame`

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot dosing (and PK) data

#### Usage

    PSSDosingData$plot(
      df_fit = NULL,
      subject_df = NULL,
      max_num_subjects = 12,
      subject_ids = NULL
    )

#### Arguments

- `df_fit`:

  Fit data frame. Uses columns `val`, `lower` and `upper`.

- `subject_df`:

  Subject data frame.

- `max_num_subjects`:

  Max number of subjects to plot.

- `subject_ids`:

  Which subjects to plot?

------------------------------------------------------------------------

### Method [`filter()`](https://rdrr.io/r/stats/filter.html)

Filter based on subject id, creates new object

#### Usage

    PSSDosingData$filter(subject_ids_keep = NULL)

#### Arguments

- `subject_ids_keep`:

  Subject ids to keep

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    PSSDosingData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
