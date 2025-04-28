# Find index of x in vector y, expecting exactly one match
find_one <- function(x, y) {
  idx <- which(y == x)
  stopifnot(length(idx) == 1)
  checkmate::assert_integerish(idx, len = 1)
  idx
}
