# Find index of x in vector y, expecting exactly one match
find_one <- function(x, y) {
  idx <- which(y == x)
  stopifnot(length(idx) == 1)
  checkmate::assert_integerish(idx, len = 1)
  idx
}

# Subject ids for simulated data
sim_subject_ids <- function(N) {
  paste0("sim-", formatC(seq_len(N), 5, flag = "0"))
}
