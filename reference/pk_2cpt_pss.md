# Partially steady-state PK model

For each subject

## Usage

``` r
pk_2cpt_pss(t, dose_ss, times, doses, theta, tau)
```

## Arguments

- t:

  vector of output time points

- dose_ss:

  dose amount in SS (for each subject)

- times:

  time points, first of which is the end of stedy-state assumption

- doses:

  doses taken after `t_last_ss`

- theta:

  PK params for each subject

- tau:

  Dosing interval (same for all subjects).

## Value

For each subject, the concentration in the central compartment at times
`t`
