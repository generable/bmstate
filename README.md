# bmstate

R package for Bayesian multistate modeling.


## Authentication

You need to set the `GITHUB_PAT` environment variable, which should be your GitHub [personal
access token](https://github.com/settings/tokens?type=beta). 

##  Installation

* Install `cmdstanr` following the instructions [here](https://mc-stan.org/cmdstanr/).
* Install `bmstate` using

```r
remotes::install_github("generable/bmstate", ref = "main", build_vignettes = TRUE)
```

Building the vignettes can take some time, so consider using `build_vignettes = FALSE`.

## Getting started

 More info is in documentation, that you can view with
```r
library(bmstate)
?bmstate
```
