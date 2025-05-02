#' Defines states and possible transitions
#'
#' @export
#' @field matrix a binary matrix where \code{matrix[i,j]} is 1 if transition
#' from state \code{i} to \code{j} is possible
#' @field states a character vector of state names
StateGraph <- R6::R6Class("StateGraph",

  # PUBLIC
  public = list(
    matrix = NULL,
    states = NULL,

    #' @description
    #' Create graph
    #'
    #' @param matrix Binary square matrix that defines the states and possible
    #' @param states Names of the states
    #' transitions
    initialize = function(matrix, states) {
      checkmate::assert_matrix(matrix, min.rows = 1, min.cols = 1)
      checkmate::assert_integerish(matrix, lower = 0, upper = 1)
      checkmate::assert_true(length(states) == nrow(matrix))
      checkmate::assert_true(nrow(matrix) == ncol(matrix))
      self$matrix <- matrix
      self$states <- states
    },

    #' Print output
    #'
    #' @return nothing
    print = function() {
      print(self$matrix)
    }
  )
)
