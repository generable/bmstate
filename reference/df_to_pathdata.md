# Create 'PathData' from a data frame of one observed path per subject

Create 'PathData' from a data frame of one observed path per subject

## Usage

``` r
df_to_pathdata(df, tm, covs = NULL, validate = TRUE)
```

## Arguments

- df:

  Data frame. Should have columns

  - `subject_id` (character)

  - `time` (numeric)

  - `state` (integer, same indexing as in `tm`)

  - `is_transition` (logical)

  - all the columns specified in `covs`

  Rules:

  - For each subject, there should be at least two rows, and the rows
    should be contiguous and the time should be non-decreasing.

  - The `is_transition` value should indicate whether the row
    corresponds to a transition.

  - The first row for each subject should never be a transition.

  - For each row that is a transition, the transition from the state of
    the previous row to the current state should be a valid transition
    in `tm`.

  - For each row that is not a transition, the state should not change
    from the previous row of that subject.

- tm:

  A
  [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)

- covs:

  covariates (character vector)

- validate:

  Do stricter data validation? Recommended to use `TRUE`.

## Value

A
[`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
object
