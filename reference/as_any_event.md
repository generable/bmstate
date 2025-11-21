# PathData to time-to-event data format for any state other than null state

PathData to time-to-event data format for any state other than null
state

## Usage

``` r
as_any_event(pd, null_state = "Randomization")
```

## Arguments

- pd:

  A
  [`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
  object

- null_state:

  Name of the base state

## Value

A
[`PathData`](https://generable.github.io/bmstate/reference/PathData.md)
object

## See also

Other PathData mutation functions:
[`as_single_event()`](https://generable.github.io/bmstate/reference/as_single_event.md),
[`as_survival()`](https://generable.github.io/bmstate/reference/as_survival.md)
