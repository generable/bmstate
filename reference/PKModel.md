# One-compartment PK Model

Class that represents a one-compartment PK model with linear oral
absorption.

## Methods

### Public methods

- [`PKModel$get_max_conc()`](#method-PKModel-get_max_conc)

- [`PKModel$set_max_conc()`](#method-PKModel-set_max_conc)

- [`PKModel$new()`](#method-PKModel-new)

- [`PKModel$ka_covs()`](#method-PKModel-ka_covs)

- [`PKModel$CL_covs()`](#method-PKModel-CL_covs)

- [`PKModel$V2_covs()`](#method-PKModel-V2_covs)

- [`PKModel$covs()`](#method-PKModel-covs)

- [`PKModel$print()`](#method-PKModel-print)

- [`PKModel$simulate_ss()`](#method-PKModel-simulate_ss)

- [`PKModel$compute_ss_auc()`](#method-PKModel-compute_ss_auc)

- [`PKModel$simulate_data()`](#method-PKModel-simulate_data)

- [`PKModel$format_params()`](#method-PKModel-format_params)

- [`PKModel$clone()`](#method-PKModel-clone)

------------------------------------------------------------------------

### Method `get_max_conc()`

Get concentration upper bound

#### Usage

    PKModel$get_max_conc()

------------------------------------------------------------------------

### Method `set_max_conc()`

Set concentration upper bound

#### Usage

    PKModel$set_max_conc(value)

#### Arguments

- `value`:

  Upper bound for concentration, to avoid numerical issues.

------------------------------------------------------------------------

### Method `new()`

Create model

#### Usage

    PKModel$new(covariates)

#### Arguments

- `covariates`:

  A list with elements `ka`, `CL`, and `V2`

------------------------------------------------------------------------

### Method `ka_covs()`

Get the covariates affecting the k_a parameter.

#### Usage

    PKModel$ka_covs()

------------------------------------------------------------------------

### Method `CL_covs()`

Get the covariates affecting the CL parameter.

#### Usage

    PKModel$CL_covs()

------------------------------------------------------------------------

### Method `V2_covs()`

Get the covariates affecting the V2 parameter.

#### Usage

    PKModel$V2_covs()

------------------------------------------------------------------------

### Method `covs()`

Get names of all unique covariates

#### Usage

    PKModel$covs()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the object info

#### Usage

    PKModel$print()

#### Returns

nothing

------------------------------------------------------------------------

### Method `simulate_ss()`

Simulate system in steady state

#### Usage

    PKModel$simulate_ss(t, theta, dose, tau)

#### Arguments

- `t`:

  output time points

- `theta`:

  parameter values

- `dose`:

  dose

- `tau`:

  dosing interval

#### Returns

A vector with the same length as `t`, representing the concentration in
the central compartment.

------------------------------------------------------------------------

### Method `compute_ss_auc()`

Compute steady-state area under curve over one dosing interval

#### Usage

    PKModel$compute_ss_auc(theta, dose)

#### Arguments

- `theta`:

  parameter values

- `dose`:

  dose

#### Returns

A numeric value

------------------------------------------------------------------------

### Method `simulate_data()`

Simulate data with many subjects

#### Usage

    PKModel$simulate_data(df_subjects, beta_pk = NULL, tau = 24, sigma = 0.3)

#### Arguments

- `df_subjects`:

  Data frame with one row for each subject

- `beta_pk`:

  Covariate effects

- `tau`:

  Dosing interval (PK)

- `sigma`:

  Noise magnitude of concentration measurements (PK)

#### Returns

Data frame with one row for each subject, and a
[`DosingData`](https://generable.github.io/bmstate/reference/DosingData.md)
object

------------------------------------------------------------------------

### Method `format_params()`

Format list of input PK parameters to standardized format.

#### Usage

    PKModel$format_params(beta_pk = NULL)

#### Arguments

- `beta_pk`:

  A list of max three elements

- `return`:

  A list with three elements

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    PKModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
