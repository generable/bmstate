#' Create a model
#'
#' @export
#' @param matrix a binary matrix where \code{matrix[i,j]} is 1 if transition
#' from state \code{i} to \code{j} is possible
#' @param states a character vector of state names
#' @param compile compile 'Stan' model?
#' @return A \code{\link{MultistateModel}} object
create_msm <- function(matrix, states, compile = TRUE) {
  tm <- TransitionMatrix$new(matrix, states)
  MultistateModel$new(tm, compile)
}

#' Main model class
#'
#' @export
#' @field transmat The transition matrix
MultistateModel <- R6::R6Class("MultistateModel",

  # PRIVATE
  private = list(
    stan_model = NULL
  ),

  # PUBLIC
  public = list(
    transmat = NULL,

    #' @description
    #' Create model
    #'
    #' @param transmat A \code{\link{TransitionMatrix}}.
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(transmat, compile = TRUE) {
      checkmate::assert_class(transmat, "TransitionMatrix")
      self$transmat <- transmat
      if (compile) {
        self$compile()
      }
    },

    #' @description Get the underlying 'Stan' model.
    get_stan_model = function() {
      private$stan_model
    },

    #' @description
    #' Create and compile the 'Stan' model.
    #'
    #' @param ... Arguments passed to \code{cmdstan_model}.
    #' @return The updated model object (invisibly).
    compile = function(...) {
      fn <- "msm.stan"
      filepath <- system.file(file.path("stan", fn), package = "bmstate")
      # silence compile warnings from cmdstan
      utils::capture.output(
        {
          mod <- cmdstanr::cmdstan_model(filepath, ...)
        },
        type = "message"
      )
      private$stan_model <- mod
      invisible(self)
    },

    #' Print the object
    #'
    #' @return nothing
    print = function() {
      self$transmat$print()
    }
  )
)
