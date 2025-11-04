# bmstate

R package for Bayesian multistate modeling.

##  Installation

* Install `cmdstanr` following the instructions [here](https://mc-stan.org/cmdstanr/).
* Install `bmstate` using

```r
remotes::install_github("generable/bmstate", ref = "main", build_vignettes = TRUE)
```

Building the vignettes can take some time, so consider using `build_vignettes = FALSE`. Use `ref = "develop"` if you want the current development version.

## Getting started

 More info is in the documentation, which you can view with
```r
library(bmstate)
?bmstate
```

See also the vignettes.

## Development

### Versioning

The package mostly conforms to [semantic versioning](https://semver.org/). The version number has the format `MAJOR.MINOR.PATCH`, where

- `MAJOR` is incremented if incompatible API changes are made
- `MINOR` version is incremented if significant new functionality is added
- `PATCH` version is incremented when making small bug fixes, very minor new functionality, improved documentation or new or improved vignettes

The version update can be done with `usethis::use_version()`.

### Development branches

When starting a new feature or a bug fix, create a branch off the `develop` branch. The version in the `develop` branch should be
something like `MAJOR.MINOR.PATCH.9000` where `MAJOR.MINOR.PATCH` is the version that exists in the `main` branch. This increment can be done with `usethis::use_dev_version()`.
When finishing the edits, create a PR to `develop`. Before the `develop` branch is merged to `main`, the patch version should be incremented. Also, make sure to run

- `styler::style_dir()` in the package root directory.
- `devtools::document()`
- `R CMD check`

### Known NOTEs in R CMD check

- The exported Stan functions seem to have no bindings in R code

```
   Undefined global functions or variables:
     compute_log_hazard_multiplier compute_theta_pk
     pop_2cpt_partly_ss_stage1 pop_2cpt_partly_ss_stage2
```

- Other

```
❯ checking package dependencies ... NOTE
  Imports includes 22 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking dependencies in R code ... NOTE
  Namespaces in Imports field not imported from:
    ‘Matrix’ ‘msm’ ‘mstate’
    All declared Imports should be used.
```

There should not be errors or warnings.

