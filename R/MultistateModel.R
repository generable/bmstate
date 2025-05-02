#' Main model class
#'
#' @export
MultiStateModel <- R6::R6Class("MultiStateModel",

  # PRIVATE
  private = list(
    stan_model = NULL
  ),

  # PUBLIC
  public = list(

    #' @description
    #' Create model
    #'
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(compile = TRUE) {
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
    #' @param dir Path to directory where to store the \code{.stan} file.
    #' @return A 'CmdStanR' model.
    create_MultiStateModel = function(dir = tempdir()) {
      code <- self$create_stancode(autoformat = FALSE)
      a <- cmdstanr::write_stan_file(code = code, dir = dir)

      # silence compile warnings from cmdstan
      utils::capture.output(
        {
          mod <- cmdstanr::cmdstan_model(a)
        },
        type = "message"
      )
      mod
    },

    #' @description
    #' Create and compile the 'Stan' model.
    #'
    #' @param ... Arguments passed to \code{create_MultiStateModel()}.
    #' @return The updated model object (invisibly).
    compile = function(...) {
      private$stan_model <- self$create_MultiStateModel(...)
      invisible(self)
    }
  )
)
