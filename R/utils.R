# Find index of x in vector y, expecting exactly one match
find_one <- function(x, y) {
  idx <- which(y == x)
  if (length(idx) != 1) {
    stop("expected exactly one match")
  }
  checkmate::assert_integerish(idx, len = 1)
  idx
}
