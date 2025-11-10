# Simulate dosing data

Simulate dosing data

## Usage

``` r
simulate_dosing(df_subjects, tau = 24, p_miss = 0.2)
```

## Arguments

- df_subjects:

  Data frame with one row for each subject

- tau:

  Dosing interval.

- p_miss:

  Probability of missing a dose.

## Value

A
[`DosingData`](https://generable.github.io/bmstate/reference/DosingData.md)
object
