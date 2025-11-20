# Path data class (R6 class)

It is not recommended for users to try to create data using the
constructor of this class. Rather use
[`df_to_pathdata`](https://generable.github.io/bmstate/reference/df_to_pathdata.md).

## Public fields

- `subject_df`:

  Data frame with one row per subject. Must have one row for each
  subject, and `subject_id` and all covariates as columns.

- `path_df`:

  Data frame of actual paths. Each row corresponds to one time point.
  Must have `path_id`, `state`, `time`, and `trans_idx` as columns.
  These are

  - `path_id`: path identifier

  - `state`: integer index of the state at the time point, and until the
    next time point

  - `time`: time point

  - `trans_idx`: Integer indicating which transition this time point
    corresponds to. If the time point is not a state transition, this
    should be 0.

- `link_df`:

  Links the path and subject data frames. Must have one row for each
  path, and `path_id`, `draw_idx`, `rep_idx`, and `subject_id` as
  columns.

- `covs`:

  Covariate column names.

- `transmat`:

  A
  [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)
  describing the system to which the paths belong.

## Methods

### Public methods

- [`PathData$new()`](#method-PathData-new)

- [`PathData$unique_subjects()`](#method-PathData-unique_subjects)

- [`PathData$covariate_names()`](#method-PathData-covariate_names)

- [`PathData$subset_covariates()`](#method-PathData-subset_covariates)

- [`PathData$state_at()`](#method-PathData-state_at)

- [`PathData$transition_times()`](#method-PathData-transition_times)

- [`PathData$n_paths()`](#method-PathData-n_paths)

- [`PathData$state_names()`](#method-PathData-state_names)

- [`PathData$get_event_state_names()`](#method-PathData-get_event_state_names)

- [`PathData$terminal_states()`](#method-PathData-terminal_states)

- [`PathData$print()`](#method-PathData-print)

- [`PathData$get_path_df()`](#method-PathData-get_path_df)

- [`PathData$as_data_frame()`](#method-PathData-as_data_frame)

- [`PathData$full_link()`](#method-PathData-full_link)

- [`PathData$as_transitions()`](#method-PathData-as_transitions)

- [`PathData$as_transitions_alt()`](#method-PathData-as_transitions_alt)

- [`PathData$as_msdata()`](#method-PathData-as_msdata)

- [`PathData$plot_paths()`](#method-PathData-plot_paths)

- [`PathData$prop_matrix()`](#method-PathData-prop_matrix)

- [`PathData$plot_graph()`](#method-PathData-plot_graph)

- [`PathData$fit_coxph()`](#method-PathData-fit_coxph)

- [`PathData$fit_mstate()`](#method-PathData-fit_mstate)

- [`PathData$filter()`](#method-PathData-filter)

- [`PathData$clone()`](#method-PathData-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize

#### Usage

    PathData$new(subject_df, path_df, link_df, transmat, covs = NULL)

#### Arguments

- `subject_df`:

  Data frame with one row per subject. Must have one row for each
  subject, and `subject_id` and all covariates as columns.

- `path_df`:

  Data frame of actual paths. Each row corresponds to one time point.
  Must have `path_id`, `state`, `time`, and `trans_idx` as columns.
  These are

  - `path_id`: path identifier

  - `state`: integer index of the state at the time point, and until the
    next time point

  - `time`: time point

  - `trans_idx`: Integer indicating which transition this time point
    corresponds to. If the time point is not a state transition, this
    should be 0.

- `link_df`:

  Links the path and subject data frames. Must have one row for each
  path, and `path_id`, `draw_idx`, `rep_idx`, and `subject_id` as
  columns.

- `transmat`:

  A
  [`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)
  describing the system to which the paths belong.

- `covs`:

  Covariate column names.

------------------------------------------------------------------------

### Method `unique_subjects()`

Get unique subject ids

#### Usage

    PathData$unique_subjects()

------------------------------------------------------------------------

### Method `covariate_names()`

Get names of covariates

#### Usage

    PathData$covariate_names()

#### Returns

a character vector

------------------------------------------------------------------------

### Method `subset_covariates()`

Create a new PathData with a subset of covariates

#### Usage

    PathData$subset_covariates(covs, renamed_old = NULL, renamed_new = NULL)

#### Arguments

- `covs`:

  Names of all new covariates

- `renamed_old`:

  Name of an old covariate to rename

- `renamed_new`:

  Name of new covariate where old one is copied

------------------------------------------------------------------------

### Method `state_at()`

For each path, get the state it is in at time t

#### Usage

    PathData$state_at(t)

#### Arguments

- `t`:

  time

#### Returns

a data frame with one row for each path

------------------------------------------------------------------------

### Method `transition_times()`

Get all transition times

#### Usage

    PathData$transition_times(trans_inds = NULL)

#### Arguments

- `trans_inds`:

  Indices of transition whose occurrence times are requested. Default is
  `NULL`, in which case all transitions are considered.

#### Returns

a numeric vector

------------------------------------------------------------------------

### Method `n_paths()`

Get number of paths

#### Usage

    PathData$n_paths()

#### Returns

an integer

------------------------------------------------------------------------

### Method `state_names()`

Get names of all states

#### Usage

    PathData$state_names()

#### Returns

character vector

------------------------------------------------------------------------

### Method `get_event_state_names()`

Get names of all states to which transitioning can be a transition event

#### Usage

    PathData$get_event_state_names()

#### Returns

character vector

------------------------------------------------------------------------

### Method `terminal_states()`

Get names of terminal states

#### Usage

    PathData$terminal_states()

#### Returns

character vector

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print info

#### Usage

    PathData$print()

#### Returns

nothing

------------------------------------------------------------------------

### Method `get_path_df()`

Get path data frame

#### Usage

    PathData$get_path_df(truncate = FALSE)

#### Arguments

- `truncate`:

  Remove rows after terminal events?

------------------------------------------------------------------------

### Method `as_data_frame()`

Convert to one long data frame

#### Usage

    PathData$as_data_frame(covariates = NULL, truncate = FALSE)

#### Arguments

- `covariates`:

  Which covariates to include?

- `truncate`:

  Remove rows after terminal events?

------------------------------------------------------------------------

### Method `full_link()`

Full link data frame

#### Usage

    PathData$full_link(covariates = NULL)

#### Arguments

- `covariates`:

  Which covariates to include

#### Returns

A `data.frame` with same number of rows as `link_df`, including also the
covariate columns and `subject_id`

------------------------------------------------------------------------

### Method `as_transitions()`

Data frame in transitions format

#### Usage

    PathData$as_transitions(covariates = NULL, truncate = FALSE)

#### Arguments

- `covariates`:

  Which covariates to include?

- `truncate`:

  Remove rows after terminal events first?

#### Returns

A `data.frame`

------------------------------------------------------------------------

### Method `as_transitions_alt()`

Data frame in alternative transitions format

#### Usage

    PathData$as_transitions_alt(covariates = NULL, truncate = FALSE)

#### Arguments

- `covariates`:

  Which covariates to include?

- `truncate`:

  Remove rows after terminal events first?

#### Returns

A `data.frame`

------------------------------------------------------------------------

### Method `as_msdata()`

Convert to format used by the 'mstate' package

#### Usage

    PathData$as_msdata(covariates = NULL)

#### Arguments

- `covariates`:

  Which covariates to include?

#### Returns

An `msdata` object

------------------------------------------------------------------------

### Method `plot_paths()`

Step plot of the paths

#### Usage

    PathData$plot_paths(n_paths = NULL, alpha = 0.5, truncate = FALSE)

#### Arguments

- `n_paths`:

  Number of paths to subsample for plotting.

- `alpha`:

  opacity

- `truncate`:

  truncate after terminal events?

------------------------------------------------------------------------

### Method `prop_matrix()`

Transition proportion matrix

#### Usage

    PathData$prop_matrix(include_censor = TRUE)

#### Arguments

- `include_censor`:

  Include censoring (no event) proportion in the matrix

#### Returns

a `table`

------------------------------------------------------------------------

### Method `plot_graph()`

Visualize the transition proportion matrix as a graph

#### Usage

    PathData$plot_graph(digits = 3, ...)

#### Arguments

- `digits`:

  Max number of digits to show in numbers

- `...`:

  Arguments passed to `qgraph`

#### Returns

`qgraph` plot

------------------------------------------------------------------------

### Method [`fit_coxph()`](https://generable.github.io/bmstate/reference/fit_coxph.md)

Fit Cox proportional hazards model

#### Usage

    PathData$fit_coxph(covariates = NULL, ...)

#### Arguments

- `covariates`:

  Covariates to include.

- `...`:

  Arguments passed to
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html).

------------------------------------------------------------------------

### Method `fit_mstate()`

Fit frequentist 'mstate' model

#### Usage

    PathData$fit_mstate(covariates = NULL, ...)

#### Arguments

- `covariates`:

  Covariates to include.

- `...`:

  Arguments passed to
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html).

------------------------------------------------------------------------

### Method [`filter()`](https://rdrr.io/r/stats/filter.html)

Filter based on subject id, creates new object

#### Usage

    PathData$filter(subject_ids_keep)

#### Arguments

- `subject_ids_keep`:

  Subject ids to keep

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    PathData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
