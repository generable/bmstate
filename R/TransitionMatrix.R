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
      checkmate::assert_true(nrow(matrix) == ncol(matrix))
      checkmate::assert_character(states, len = nrow(matrix), unique = TRUE)
      self$matrix <- matrix
      self$states <- states
    },

    #' Number of states
    #'
    #' @return integer \code{N}
    num_states = function() {
      nrow(self$matrix)
    },

    #' Number of transitions
    #'
    #' @return integer
    num_trans = function() {
      sum(as.vector(self$matrix))
    },

    #' @description Get number of different transition types.
    num_trans_types = function() {
      length(unique(self$trans_df()$trans_type))
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

    #' @description Transition index matrix to format used by 'mstate'
    #'
    #' @return a matrix
    as_mstate_transmat = function() {
      TFI <- self$as_transition_index_matrix()
      TFI[TFI == 0] <- NA
      TFI
    },

    #' @description Get indices of possible transitions from given state
    #' @param state Index of state
    possible_transitions_from = function(state) {
      S <- self$num_states()
      checkmate::assert_integerish(state, len = 1, lower = 1, upper = S)
      row <- self$as_transition_index_matrix()[state, ]
      as.numeric(row[row > 0])
    },

    #' Get target state of a transition
    #'
    #' @param trans_idx Index of the transition
    target_state = function(trans_idx) {
      H <- self$num_trans()
      checkmate::assert_integerish(trans_idx, len = 1, lower = 1, upper = H)
      self$trans_df()$state[trans_idx]
    },

    #' Get states that are at risk when at given state
    #'
    #' @param state index of current state
    #' @return indices of states at risk
    at_risk = function(state) {
      possible_trans <- self$possible_transitions_from(state)
      if (length(possible_trans) == 0) {
        return(integer(0))
      }
      self$trans_df()[possible_trans, ]$state
    },

    #' @description A data frame of states
    #'
    #' @return a \code{data.frame}
    states_df = function() {
      N <- self$num_states()
      idx <- seq_len(N)
      abs <- rep(FALSE, N)
      src <- rep(FALSE, N)
      abs[self$absorbing_states(names = FALSE)] <- TRUE
      src[self$source_states(names = FALSE)] <- TRUE
      data.frame(
        state_idx = idx, state = self$states, terminal = abs, source = src
      )
    },

    #' @description A data frame of transitions ("legend")
    #'
    #' @return a \code{data.frame}
    trans_df = function() {
      H <- self$num_trans()
      TFI <- self$as_transition_index_matrix()
      idx <- seq_len(H)
      prev_state <- rep(0, H)
      target_state <- rep(0, H)
      for (i in 1:nrow(TFI)) {
        for (j in 1:ncol(TFI)) {
          if (TFI[i, j] > 0) {
            prev_state[TFI[i, j]] <- i
            target_state[TFI[i, j]] <- j
          }
        }
      }
      df <- data.frame(
        trans_idx = idx,
        prev_state = prev_state,
        state = target_state
      )
      df$trans_char <- paste0(
        self$states[df$prev_state], " -> ",
        self$states[df$state]
      )
      df$trans_type <- dplyr::dense_rank(df$state)
      df
    },

    #' @description Visualize the matrix as a graph
    #'
    #' @param edge_labs Edge labels
    #' @param ... Arguments passed to \code{qgraph::qgraph}
    #' @return \code{qgraph} plot
    plot = function(edge_labs = FALSE, ...) {
      transition_matrix_plot(
        self$as_matrix(),
        self$absorbing_states(),
        self$source_states(),
        edge_labs = edge_labs,
        ...
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
      out <- which(colSums(self$matrix) == 0)
      if (names) {
        out <- self$states[out]
      }
      out
    }
  )
)


#' Create a common transition matrix
#'
#' @description
#' \itemize{
#' \item \code{full}: All possible transitions
#' between non-terminal and non-source states are added, as well as transition
#' to terminal states from all states.
#' \item \code{survival}: Basic survival model with single transition
#' \item \code{comprisk}: Competing risks transition matrix
#' \item \code{illnessdeath}: Illness-Death model
#' \item \code{progression}: Disease progression model
#' \item \code{diamond}: Two intermediate states both leading to same terminal
#' state
#' }
#' @param state_names Names of the states. The length of this character
#' vector defines the number of states.
#' @param sources Indices of source states.
#' @param terminal Indices of terminal states
#' @param self_loops Add self-loops to non-terminal and non-source states?
#' @return A \code{\link{TransitionMatrix}}
#' @name transmat
NULL

#' @export
#' @rdname transmat
transmat_full <- function(state_names = LETTERS[1:4], sources = 1,
                          terminal = length(state_names), self_loops = TRUE) {
  N <- length(state_names)
  checkmate::assert_character(state_names, min.len = 1)
  checkmate::assert_integerish(sources, lower = 1, upper = N)
  checkmate::assert_integerish(terminal, lower = 1, upper = N)
  checkmate::assert_logical(self_loops, len = 1)
  mat <- matrix(1, N, N)
  mat[, sources] <- 0
  mat[terminal, ] <- 0
  if (!self_loops) {
    for (n in seq_len(N)) {
      mat[n, n] <- 0
    }
  }
  TransitionMatrix$new(mat, state_names)
}

#' @export
#' @rdname transmat
transmat_comprisk <- function(state_names = LETTERS[1:4]) {
  N <- length(state_names)
  checkmate::assert_character(state_names, min.len = 3)
  mat <- matrix(0, N, N)
  mat[1, 2:(ncol(mat))] <- 1
  TransitionMatrix$new(mat, state_names)
}

#' @export
#' @rdname transmat
transmat_survival <- function(state_names = c("Alive", "Dead")) {
  checkmate::assert_character(state_names, len = 2)
  transmat_full(state_names, self_loops = F)
}

#' @export
#' @rdname transmat
transmat_illnessdeath <- function(state_names =
                                    c("Healthy", "Diseased", "Dead")) {
  checkmate::assert_character(state_names, len = 3)
  transmat_full(state_names, self_loops = F)
}

#' @export
#' @rdname transmat
transmat_progression <- function(state_names =
                                   c("Healthy", "Mild", "Severe", "Dead")) {
  mat <- matrix(0, 4, 4)
  mat[1, 2] <- 1 # Healthy to Mild
  mat[1, 3] <- 1 # Healthy to Severe
  mat[2, 3] <- 1 # Mild to Severe
  mat[c(1, 2, 3), 4] <- 1 # Any -> Dead
  TransitionMatrix$new(mat, state_names)
}

#' @export
#' @rdname transmat
transmat_diamond <- function(state_names = LETTERS[1:4]) {
  checkmate::assert_character(state_names, len = 4)
  mat <- matrix(0, 4, 4)
  mat[1, 2] <- 1
  mat[2, 3] <- 1
  mat[3, 2] <- 1
  mat[3, 4] <- 1
  mat[2, 4] <- 1
  mat[1, 4] <- 1
  mat[1, 3] <- 1
  mat[2, 2] <- 0
  mat[3, 3] <- 0
  TransitionMatrix$new(mat, state_names)
}


# Check that two transition matrices are the same
check_equal_transmats <- function(tm1, tm2) {
  stopifnot(all(tm1$states == tm2$states))
  stopifnot(isTRUE(all.equal(tm1$matrix, tm2$matrix)))
}
