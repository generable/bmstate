# Find index of x in vector y, expecting exactly one match
find_one <- function(x, y) {
  idx <- which(y == x)
  stopifnot(length(idx) == 1)
  checkmate::assert_integerish(idx, len = 1)
  idx
}


#' Train-test split
#'
#' @export
#' @param pd A \code{\link{PathData}} object
#' @param p_test Proportion of test subjects
#' @return a list with test and train subject indices
do_split <- function(pd, p_test = 0.25) {
  checkmate::assert_class(pd, "PathData")
  checkmate::assert_number(p_test, lower = 0, upper = 1)
  all_sub <- unique(as.numeric(pd$subject_df$subject_index))
  N_sub <- length(all_sub)
  i_test <- sample.int(N_sub, round(N_sub * p_test))
  test_sub <- all_sub[i_test]
  train_sub <- setdiff(all_sub, test_sub)
  list(
    test_sub = test_sub,
    train_sub = train_sub
  )
}
