#' Defines states and possible transitions
#'
#' @export
#' @field matrix a binary \code{N} x \code{N} matrix where
#' \code{matrix[i,j]} is 1 if transition
#' from state \code{i} to \code{j} is possible
#' @field states a character vector of state names (length \code{N})
TransitionMatrix <- R6::R6Class("TransitionMatrix",

  # PUBLIC
  public = list(
    matrix = NULL,
    states = NULL,

    #' @description
    #' Create graph
    #'
    #' @param matrix a binary \code{N} x \code{N} matrix where
    #' \code{matrix[i,j]} is 1 if transition
    #' from state \code{i} to \code{j} is possible
    #' @param states a character vector of state names (length \code{N})
    initialize = function(matrix, states) {
      checkmate::assert_matrix(matrix, min.rows = 1, min.cols = 1)
      checkmate::assert_integerish(matrix, lower = 0, upper = 1)
      checkmate::assert_true(length(states) == nrow(matrix))
      checkmate::assert_true(nrow(matrix) == ncol(matrix))
      self$matrix <- matrix
      self$states <- states
    },

    #' Number of states
    #'
    #' @return integer \code{N}
    num_states = function() {
      nrow(self$matrix)
    },

    #' @description As a matrix
    #'
    #' @return An \code{N} x \code{N} matrix with named columns and rows
    as_matrix = function() {
      out <- self$matrix
      colnames(out) <- self$states
      rownames(out) <- colnames(out)
      out
    },

    #' @description As a matrix indexing the transitions
    #'
    #' @return An \code{N} x \code{N} matrix with named columns and rows
    as_transition_index_matrix = function() {
      a <- self$as_matrix()
      out <- matrix(0, nrow(a), ncol(a))
      count <- 0
      for (i in 1:nrow(a)) {
        for (j in 1:ncol(a)) {
          if (a[i, j] == 1) {
            count <- count + 1
            out[i, j] <- count
          }
        }
      }
      colnames(out) <- colnames(a)
      rownames(out) <- rownames(a)
      out
    },

    #' @description A data frame of states
    #'
    #' @return a \code{data.frame} with columns \code{state_idx} and \code{state}
    states_df = function() {
      N <- self$num_states()
      idx <- seq_len(N)
      abs <- rep(FALSE, N)
      src <- rep(FALSE, N)
      abs[self$absorbing_states(names = FALSE)] <- TRUE
      src[self$source_states(names = FALSE)] <- TRUE
      data.frame(
        state_idx = idx, state = self$states, absorbing = abs,
        source = src
      )
    },

    #' @description Visualize the matrix as a graph
    #'
    #' @param ... Arguments passed to \code{qgraph}
    #' @return \code{qgraph} plot
    plot = function(...) {
      transition_matrix_plot(
        self$as_matrix(),
        self$absorbing_states(),
        self$source_states(),
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
