# Simulate dosing data

Simulate dosing data

## Usage

``` r
simulate_dosing(df_subjects, tau = 24, p_miss = 0.2, t_jitter = 4)
```

## Arguments

- df_subjects:

  Data frame with one row for each subject. Must have columns
  `subject_id, num_ss_doses, num_doses, dose`.

- tau:

  Supposed dosing interval (same for each subject).

- p_miss:

  Probability of missing a dose.

- t_jitter:

  Randomness added to dose times.

## Value

A
[`DosingData`](https://generable.github.io/bmstate/reference/DosingData.md)
object
