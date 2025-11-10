# Joint data class (R6 class)

Joint data class (R6 class)

Joint data class (R6 class)

## Public fields

- `paths`:

  The events data, an object of class
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md).

- `dosing`:

  The dosing data, an object of class
  [`DosingData`](https://generable.github.io/bmstate/reference/DosingData.md).
  Initialize

## Methods

### Public methods

- [`JointData$new()`](#method-JointData-new)

- [`JointData$filter()`](#method-JointData-filter)

- [`JointData$print()`](#method-JointData-print)

- [`JointData$plot_dosing()`](#method-JointData-plot_dosing)

- [`JointData$clone()`](#method-JointData-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    JointData$new(paths, dosing)

#### Arguments

- `paths`:

  The events data, an object of class
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md).

- `dosing`:

  The dosing data, an object of class
  [`DosingData`](https://generable.github.io/bmstate/reference/DosingData.md).
  Can be also `NULL`. Filter subjects, creates new object

------------------------------------------------------------------------

### Method [`filter()`](https://rdrr.io/r/stats/filter.html)

#### Usage

    JointData$filter(subject_ids_keep)

#### Arguments

- `subject_ids_keep`:

  Subjects to keep

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print info

#### Usage

    JointData$print()

------------------------------------------------------------------------

### Method `plot_dosing()`

Plot dosing

#### Usage

    JointData$plot_dosing(df_fit = NULL, max_num_subjects = 12)

#### Arguments

- `df_fit`:

  Fit data frame.

- `max_num_subjects`:

  Max number of subjects to plot.

#### Returns

a `ggplot`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    JointData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
