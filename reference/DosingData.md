# Dosing data class (R6 class)

Dosing data class (R6 class)

Dosing data class (R6 class)

## Public fields

- `doses`:

  Dose amounts (list with length equal to number of subjects).

- `times`:

  Dose times (list with length equal to number of subjects).

- `dose_ss`:

  Steady-state dose.

- `tau_ss`:

  Steady-state dosing interval.

- `subject_ids`:

  Subject ids. Initialize

## Methods

### Public methods

- [`DosingData$new()`](#method-DosingData-new)

- [`DosingData$as_data_frame()`](#method-DosingData-as_data_frame)

- [`DosingData$print()`](#method-DosingData-print)

- [`DosingData$num_subjects()`](#method-DosingData-num_subjects)

- [`DosingData$simulate_pk()`](#method-DosingData-simulate_pk)

- [`DosingData$plot()`](#method-DosingData-plot)

- [`DosingData$filter()`](#method-DosingData-filter)

- [`DosingData$clone()`](#method-DosingData-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    DosingData$new(subject_ids, doses, times, dose_ss = NULL, tau_ss = 24)

#### Arguments

- `subject_ids`:

  A character vector

- `doses`:

  Dose amounts (list with length equal to number of subjects).

- `times`:

  Dose times (list with length equal to number of subjects).

- `dose_ss`:

  Steady-state dose.

- `tau_ss`:

  Steady-state dosing interval.

------------------------------------------------------------------------

### Method `as_data_frame()`

As data frame

#### Usage

    DosingData$as_data_frame()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print info

#### Usage

    DosingData$print()

------------------------------------------------------------------------

### Method `num_subjects()`

Get number of subjects Simulate PK dynamics

#### Usage

    DosingData$num_subjects()

------------------------------------------------------------------------

### Method `simulate_pk()`

#### Usage

    DosingData$simulate_pk(t, theta)

#### Arguments

- `t`:

  A vector of output times for each subject (a list).

- `theta`:

  A matrix of parameters.

#### Returns

a `data.frame` Plot dosing (and PK) data

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

#### Usage

    DosingData$plot(df_fit = NULL, subject_df = NULL, max_num_subjects = 12)

#### Arguments

- `df_fit`:

  Fit data frame.

- `subject_df`:

  Subject data frame.

- `max_num_subjects`:

  Max number of subjects to plot.

------------------------------------------------------------------------

### Method [`filter()`](https://rdrr.io/r/stats/filter.html)

Filter based on subject id, creates new object

#### Usage

    DosingData$filter(subject_ids_keep = NULL)

#### Arguments

- `subject_ids_keep`:

  Subject ids to keep

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DosingData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
