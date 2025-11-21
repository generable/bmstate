# PathData to time-to-event data format with a single event

PathData to time-to-event data format with a single event

## Usage

``` r
as_single_event(pd, event, null_state = "Randomization")
```

## Arguments

- pd:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object

- event:

  Name of the state corresponding to the event of interest (character)

- null_state:

  Name of the base state

## Value

A
[`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
object

## See also

Other PathData mutation functions:
[`as_any_event()`](https://generable.github.io/bmstate/reference/as_any_event.md),
[`as_survival()`](https://generable.github.io/bmstate/reference/as_survival.md)
