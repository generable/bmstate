# Create a data frame to evaluate state visit probability predictions

Create a data frame to evaluate state visit probability predictions

## Usage

``` r
create_scoring_df(paths_obs, paths_pred, state_name)
```

## Arguments

- paths_obs:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object of observed paths.

- paths_pred:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object of predicted paths.

- state_name:

  Name of the state of interest.

## Value

A `data.frame` with `surv` and `pred_prob` columns
