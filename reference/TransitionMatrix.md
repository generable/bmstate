# Defines states and possible transitions

Defines states and possible transitions

Defines states and possible transitions

## Public fields

- `matrix`:

  a binary `N` x `N` matrix where `matrix[i,j]` is 1 if transition from
  state `i` to `j` is possible

- `states`:

  a character vector of state names (length `N`)

## Methods

### Public methods

- [`TransitionMatrix$new()`](#method-TransitionMatrix-new)

- [`TransitionMatrix$num_states()`](#method-TransitionMatrix-num_states)

- [`TransitionMatrix$num_trans()`](#method-TransitionMatrix-num_trans)

- [`TransitionMatrix$num_trans_types()`](#method-TransitionMatrix-num_trans_types)

- [`TransitionMatrix$as_matrix()`](#method-TransitionMatrix-as_matrix)

- [`TransitionMatrix$as_transition_index_matrix()`](#method-TransitionMatrix-as_transition_index_matrix)

- [`TransitionMatrix$as_mstate_transmat()`](#method-TransitionMatrix-as_mstate_transmat)

- [`TransitionMatrix$possible_transitions_from()`](#method-TransitionMatrix-possible_transitions_from)

- [`TransitionMatrix$target_state()`](#method-TransitionMatrix-target_state)

- [`TransitionMatrix$at_risk()`](#method-TransitionMatrix-at_risk)

- [`TransitionMatrix$states_df()`](#method-TransitionMatrix-states_df)

- [`TransitionMatrix$trans_df()`](#method-TransitionMatrix-trans_df)

- [`TransitionMatrix$plot()`](#method-TransitionMatrix-plot)

- [`TransitionMatrix$print()`](#method-TransitionMatrix-print)

- [`TransitionMatrix$absorbing_states()`](#method-TransitionMatrix-absorbing_states)

- [`TransitionMatrix$source_states()`](#method-TransitionMatrix-source_states)

- [`TransitionMatrix$clone()`](#method-TransitionMatrix-clone)

------------------------------------------------------------------------

### Method `new()`

Create graph

#### Usage

    TransitionMatrix$new(matrix, states)

#### Arguments

- `matrix`:

  a binary `N` x `N` matrix where `matrix[i,j]` is 1 if transition from
  state `i` to `j` is possible

- `states`:

  a character vector of state names (length `N`) Number of states

------------------------------------------------------------------------

### Method `num_states()`

#### Usage

    TransitionMatrix$num_states()

#### Returns

integer `N` Number of transitions

------------------------------------------------------------------------

### Method `num_trans()`

#### Usage

    TransitionMatrix$num_trans()

#### Returns

integer

------------------------------------------------------------------------

### Method `num_trans_types()`

Get number of different transition types.

#### Usage

    TransitionMatrix$num_trans_types()

------------------------------------------------------------------------

### Method `as_matrix()`

As a matrix

#### Usage

    TransitionMatrix$as_matrix()

#### Returns

An `N` x `N` matrix with named columns and rows

------------------------------------------------------------------------

### Method `as_transition_index_matrix()`

As a matrix indexing the transitions

#### Usage

    TransitionMatrix$as_transition_index_matrix()

#### Returns

An `N` x `N` matrix with named columns and rows

------------------------------------------------------------------------

### Method `as_mstate_transmat()`

Transition index matrix to format used by 'mstate'

#### Usage

    TransitionMatrix$as_mstate_transmat()

#### Returns

a matrix

------------------------------------------------------------------------

### Method `possible_transitions_from()`

Get indices of possible transitions from given state

#### Usage

    TransitionMatrix$possible_transitions_from(state)

#### Arguments

- `state`:

  Index of state Get target state of a transition

------------------------------------------------------------------------

### Method `target_state()`

#### Usage

    TransitionMatrix$target_state(trans_idx)

#### Arguments

- `trans_idx`:

  Index of the transition Get states that are at risk when at given
  state

------------------------------------------------------------------------

### Method `at_risk()`

#### Usage

    TransitionMatrix$at_risk(state)

#### Arguments

- `state`:

  index of current state

#### Returns

indices of states at risk

------------------------------------------------------------------------

### Method `states_df()`

A data frame of states

#### Usage

    TransitionMatrix$states_df()

#### Returns

a `data.frame`

------------------------------------------------------------------------

### Method `trans_df()`

A data frame of transitions ("legend")

#### Usage

    TransitionMatrix$trans_df()

#### Returns

a `data.frame`

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Visualize the matrix as a graph

#### Usage

    TransitionMatrix$plot(edge_labs = FALSE, ...)

#### Arguments

- `edge_labs`:

  Edge labels

- `...`:

  Arguments passed to
  [`qgraph::qgraph`](https://rdrr.io/pkg/qgraph/man/qgraph.html)

#### Returns

`qgraph` plot Print output

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    TransitionMatrix$print()

#### Returns

nothing Get states that cannot be transitioned from

------------------------------------------------------------------------

### Method `absorbing_states()`

#### Usage

    TransitionMatrix$absorbing_states(names = TRUE)

#### Arguments

- `names`:

  Return names of the states? Otherwise returns indices. Get states that
  cannot be transitioned to

------------------------------------------------------------------------

### Method `source_states()`

#### Usage

    TransitionMatrix$source_states(names = TRUE)

#### Arguments

- `names`:

  Return names of the states? Otherwise returns indices.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    TransitionMatrix$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
