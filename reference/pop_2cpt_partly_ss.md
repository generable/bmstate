# Partially steady-state PK model

For each subject

## Usage

``` r
pop_2cpt_partly_ss(t, dose_ss, times, doses, theta, tau)
```

## Arguments

- t:

  vector of output time points

- dose_ss:

  dose amount in SS (for each subject)

- times:

  time points after `t_last_ss`

- doses:

  doses taken after `t_last_ss`

- theta:

  PK params for each subject

- tau:

  Dosing interval (same for all subjects).
