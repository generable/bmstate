# For each non-source state, compute probability of visiting it at least once before given time

For each non-source state, compute probability of visiting it at least
once before given time

## Usage

``` r
p_state_visit(pd, t = NULL, by = NULL)
```

## Arguments

- pd:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object.

- t:

  The given time. If `NULL`, is set to `max(pd$get_path_df()$time)`.

- by:

  Factor to summarize over.

## Value

A data frame

## See also

Other PathData summary functions:
[`p_state_visit_per_subject()`](https://generable.github.io/bmstate/reference/p_state_visit_per_subject.md)
