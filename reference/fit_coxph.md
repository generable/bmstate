# Fit Cox PH model using Breslow method

Fit Cox PH model using Breslow method

## Usage

``` r
fit_coxph(msdat, formula_rhs = NULL)
```

## Arguments

- msdat:

  `msdata` object

- formula_rhs:

  Formula right hand side that is appended to
  `Surv(Tstart, Tstop, status) ~ `. If `NULL` (default), then
  `strata(trans)` is used

## Value

value returned by
[`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
