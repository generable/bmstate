# Create a common transition matrix

- `full`: All possible transitions between non-terminal and non-source
  states are added, as well as transition to terminal states from all
  states.

- `survival`: Basic survival model with single transition

- `comprisk`: Competing risks transition matrix

- `illnessdeath`: Illness-Death model

- `progression`: Disease progression model

- `diamond`: Two intermediate states both leading to same terminal state

## Usage

``` r
transmat_full(
  state_names = LETTERS[1:4],
  sources = 1,
  terminal = length(state_names),
  self_loops = TRUE
)

transmat_comprisk(state_names = LETTERS[1:4])

transmat_survival(state_names = c("Alive", "Dead"))

transmat_illnessdeath(state_names = c("Healthy", "Diseased", "Dead"))

transmat_progression(state_names = c("Healthy", "Mild", "Severe", "Dead"))

transmat_diamond(state_names = LETTERS[1:4])
```

## Arguments

- state_names:

  Names of the states. The length of this character vector defines the
  number of states.

- sources:

  Indices of source states.

- terminal:

  Indices of terminal states

- self_loops:

  Add self-loops to non-terminal and non-source states?

## Value

A
[`TransitionMatrix`](https://generable.github.io/bmstate/reference/TransitionMatrix.md)
