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

When starting a new feature or a bug fix, create a branch off the `develop` branch. The version in the `development` branch should be
something like `MAJOR.MINOR.PATCH.9000` where `MAJOR.MINOR.PATCH` is the version that exists in the `main` branch. This increment can be done with `usethis::use_dev_version()`.
When finishing the edits, create a PR to `develop`. Before the `develop` branch is merged to `main`, the patch version should be incremented. Also, make sure to run

- `styler::style_dir()` in the package root directory.
- `devtools::document()`
- `R CMD check`

### Known issues in R CMD check

- A warning can be raised because of a too small `delta_grid` which is due to randomness in simulation of test data for unit tests
- The exported Stan functions seem to have no bindings in R code

These are acceptable warnings/notes that can be ignored for now.
