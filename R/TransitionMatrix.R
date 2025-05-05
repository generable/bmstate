#' Defines states and possible transitions
#'
#' @export
#' @field matrix a binary \code{N} x \code{N} matrix where
#' \code{matrix[i,j]} is 1 if transition
#' from state \code{i} to \code{j} is possible
#' @field states a character vector of state names (length \code{N})
#' @field censor_state name of a censoring state (not included in the matrix
#' or the state names)
TransitionMatrix <- R6::R6Class("TransitionMatrix",

  # PUBLIC
  public = list(
    matrix = NULL,
    states = NULL,
    censor_state = NULL,

    #' @description
    #' Create graph
    #'
    #' @param matrix a binary \code{N} x \code{N} matrix where
    #' \code{matrix[i,j]} is 1 if transition
    #' from state \code{i} to \code{j} is possible
    #' @param states a character vector of state names (length \code{N})
    #' @param censor_state name of a censoring state (not included in the matrix
    #' or the state names)
    #' transitions
    initialize = function(matrix, states, censor_state = "Censor") {
      checkmate::assert_matrix(matrix, min.rows = 1, min.cols = 1)
      checkmate::assert_integerish(matrix, lower = 0, upper = 1)
      checkmate::assert_true(length(states) == nrow(matrix))
      checkmate::assert_true(nrow(matrix) == ncol(matrix))
      self$matrix <- matrix
      self$states <- states
      self$censor_state <- censor_state
    },

    #' Number of states, not including censor
    #'
    #' @return integer \code{N}
    num_states = function() {
      nrow(self$matrix)
    },

    #' @description As a matrix
    #'
    #' @return An \code{N+1} x \code{N+1} matrix with named columns and rows
    as_matrix = function() {
      r <- self$matrix
      N <- self$num_states()
      out <- matrix(0, N + 1, N + 1)
      out[1:N, 1:N] <- r
      out[, N + 1] <- 1 # edge from every state to censor
      out[N + 1, N + 1] <- 0 # no loop from censor to itself
      colnames(out) <- c(self$states, self$censor_state)
      rownames(out) <- colnames(out)
      out
    },

    #' @description Visualize the matrix as a graph
    #'
    #' @param include_censor Include censoring state in the graph?
    #' @param ... Arguments passed to \code{qgraph}
    #' @return \code{qgraph} plot
    plot = function(include_censor = FALSE, ...) {
      f <- self$as_matrix()
      transition_matrix_plot(
        f, self$absorbing_states(),
        self$censor_state, self$source_states(),
        include_censor,
        edge_labs = FALSE, ...
      )
    },

    #' Print output
    #'
    #' @return nothing
    print = function() {
      print(self$as_matrix())
    },

    #' Get states that cannot be transitioned from
    #'
    #' @param names Return names of the states? Otherwise returns indices.
    absorbing_states = function(names = TRUE) {
      out <- which(rowSums(self$matrix) == 0)
      if (names) {
        out <- self$states[out]
      }
      out
    },

    #' Get states that cannot be transitioned to
    #'
    #' @param names Return names of the states? Otherwise returns indices.
    source_states = function(names = TRUE) {
      out <- which(colSums(s$matrix) == 0)
      if (names) {
        out <- self$states[out]
      }
      out
    }
  )
)
