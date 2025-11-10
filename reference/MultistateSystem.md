# A multistate system with proportional spline hazards

A multistate system with proportional spline hazards

A multistate system with proportional spline hazards

## Methods

### Public methods

- [`MultistateSystem$new()`](#method-MultistateSystem-new)

- [`MultistateSystem$num_states()`](#method-MultistateSystem-num_states)

- [`MultistateSystem$num_trans()`](#method-MultistateSystem-num_trans)

- [`MultistateSystem$num_weights()`](#method-MultistateSystem-num_weights)

- [`MultistateSystem$num_knots()`](#method-MultistateSystem-num_knots)

- [`MultistateSystem$tm()`](#method-MultistateSystem-tm)

- [`MultistateSystem$knots_not_set()`](#method-MultistateSystem-knots_not_set)

- [`MultistateSystem$set_knots()`](#method-MultistateSystem-set_knots)

- [`MultistateSystem$get_knots()`](#method-MultistateSystem-get_knots)

- [`MultistateSystem$get_tmax()`](#method-MultistateSystem-get_tmax)

- [`MultistateSystem$print()`](#method-MultistateSystem-print)

- [`MultistateSystem$log_baseline_hazard()`](#method-MultistateSystem-log_baseline_hazard)

- [`MultistateSystem$basisfun_matrix()`](#method-MultistateSystem-basisfun_matrix)

- [`MultistateSystem$log_inst_hazard()`](#method-MultistateSystem-log_inst_hazard)

- [`MultistateSystem$intensity_matrix()`](#method-MultistateSystem-intensity_matrix)

- [`MultistateSystem$has_self_loops()`](#method-MultistateSystem-has_self_loops)

- [`MultistateSystem$max_inst_hazard()`](#method-MultistateSystem-max_inst_hazard)

- [`MultistateSystem$simulate()`](#method-MultistateSystem-simulate)

- [`MultistateSystem$clone()`](#method-MultistateSystem-clone)

------------------------------------------------------------------------

### Method `new()`

Create model

#### Usage

    MultistateSystem$new(transmat)

#### Arguments

- `transmat`:

  A
  [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md).

------------------------------------------------------------------------

### Method `num_states()`

Get number of states

#### Usage

    MultistateSystem$num_states()

#### Returns

integer

------------------------------------------------------------------------

### Method `num_trans()`

Get number of transitions

#### Usage

    MultistateSystem$num_trans()

#### Returns

integer Get number of spline weight parameters

------------------------------------------------------------------------

### Method `num_weights()`

#### Usage

    MultistateSystem$num_weights()

#### Returns

an integer (`num_knots` - 2) + `self$spline_k`

------------------------------------------------------------------------

### Method `num_knots()`

Get number of spline knots

#### Usage

    MultistateSystem$num_knots()

#### Returns

integer, zero if knots have not been set

------------------------------------------------------------------------

### Method `tm()`

Get transition matrix

#### Usage

    MultistateSystem$tm()

------------------------------------------------------------------------

### Method `knots_not_set()`

Have knots not been set yet?

#### Usage

    MultistateSystem$knots_not_set()

#### Returns

a logical

------------------------------------------------------------------------

### Method `set_knots()`

Set spline knot locations

#### Usage

    MultistateSystem$set_knots(locations)

#### Arguments

- `locations`:

  A numeric vector

------------------------------------------------------------------------

### Method `get_knots()`

Get knot locations

#### Usage

    MultistateSystem$get_knots()

#### Returns

a numeric vector

------------------------------------------------------------------------

### Method `get_tmax()`

Get max time set for the model

#### Usage

    MultistateSystem$get_tmax()

#### Returns

a number Print the object

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    MultistateSystem$print()

#### Returns

nothing

------------------------------------------------------------------------

### Method `log_baseline_hazard()`

Evaluate log baseline hazard

#### Usage

    MultistateSystem$log_baseline_hazard(t, log_w0, w = NULL, SBF = NULL)

#### Arguments

- `t`:

  Output time points. Not used if `SBF` is given.

- `log_w0`:

  Intercept (log_scale).

- `w`:

  Spline weights (log scale). If `NULL`, will be set to a vector of
  zeros, meaning that the log hazard is constant at `log_w0`

- `SBF`:

  Precomputed basisfunction matrix at `t`.

#### Returns

a vector with same length as `t`

------------------------------------------------------------------------

### Method `basisfun_matrix()`

Evaluate basis function matrix for baseline hazard

#### Usage

    MultistateSystem$basisfun_matrix(t)

#### Arguments

- `t`:

  Evaluation time points.

#### Returns

Matrix with shape `c(N, L+1)` where `L` is number of knots

------------------------------------------------------------------------

### Method `log_inst_hazard()`

Evaluate log instant hazard

#### Usage

    MultistateSystem$log_inst_hazard(t, w, log_w0, log_m, SBF = NULL)

#### Arguments

- `t`:

  Time point(s). Not used if `SBF` is given.

- `w`:

  Spline basis function weights (vector)

- `log_w0`:

  Intercept (log)

- `log_m`:

  Hazard multiplier (log)

- `SBF`:

  Pre-computed basis function matrix at `t`.

------------------------------------------------------------------------

### Method `intensity_matrix()`

Evaluate transition intensity (generator) matrix at time t

#### Usage

    MultistateSystem$intensity_matrix(t, log_w0, w = NULL, log_m = NULL)

#### Arguments

- `t`:

  A number (time point)

- `log_w0`:

  A vector of length `n_trans`

- `w`:

  An array of shape `n_trans` x `n_weights`

- `log_m`:

  A vector of length `n_trans`

#### Returns

a matrix with shape `n_states` x `n_states`

------------------------------------------------------------------------

### Method `has_self_loops()`

Does the system have self loops?

#### Usage

    MultistateSystem$has_self_loops()

#### Returns

Boolean value

------------------------------------------------------------------------

### Method `max_inst_hazard()`

Max instant hazard on interval (t1, t2)

#### Usage

    MultistateSystem$max_inst_hazard(t1, t2, w, log_w0, log_m)

#### Arguments

- `t1`:

  Start time point

- `t2`:

  End time point

- `w`:

  Spline basis function weights (vector)

- `log_w0`:

  Intercept (log)

- `log_m`:

  Hazard multiplier (log) Generate paths

------------------------------------------------------------------------

### Method [`simulate()`](https://rdrr.io/r/stats/simulate.html)

#### Usage

    MultistateSystem$simulate(
      w,
      log_w0,
      log_m,
      init_state = 1,
      t_start = 0,
      t_max = NULL,
      n_rep = 1
    )

#### Arguments

- `w`:

  An array of shape `n_draws` x `n_trans` x `n_weights`

- `log_w0`:

  An array of shape `n_draws` x `n_trans`

- `log_m`:

  An array of shape `n_draws` x `n_trans`

- `init_state`:

  Index of starting state. A single value or a vector with length equal
  to `n_draws`.

- `t_start`:

  Start time.

- `t_max`:

  Max time. If `NULL`, the max time of the model is used.

- `n_rep`:

  Number of repetitions to do for each draw.

#### Returns

A data frame with `n_draws` x `n_rep` paths.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultistateSystem$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
