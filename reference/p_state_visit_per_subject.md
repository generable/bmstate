# For each subject, compute probability of visiting a given state at least once before given time

Convenient wrapper for
[`p_state_visit`](https://generable.github.io/bmstate/reference/p_state_visit.md).

## Usage

``` r
p_state_visit_per_subject(pd, state_name, t = NULL)
```

## Arguments

- pd:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object.

- state_name:

  Name of the state (character).

- t:

  The given time. If `NULL`, is set to `max(pd$get_path_df()$time)`.

## Value

A data frame

## See also

Other PathData summary functions:
[`p_state_visit()`](https://generable.github.io/bmstate/reference/p_state_visit.md)
