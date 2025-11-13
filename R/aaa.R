#' The 'bmstate' package.
#'
#' @description Bayesian multistate modeling of right-censored time-to-event data.
#' @author Juho Timonen and Jacqueline Buros
#' @keywords multistate Bayesian
#'
#' @section Transition matrix:
#' For any analysis with this package, you will need a \code{\link{TransitionMatrix}}
#' which specifies the states and possible transitions between them.
#' This can be created directly with
#' \code{TransitionMatrix$new()} or using built-in functions for common
#' multistate systems such as \code{\link{transmat_illnessdeath}}.
#'
#' @section Data format:
#' See the function \code{\link{df_to_pathdata}} for how to create a
#' \code{\link{PathData}} object from your data frame and a given
#' \code{\link{TransitionMatrix}}.
#'
#' @section Creating a model:
#' You can create a model using
#' \code{\link{create_msm}}.
#' A model is represented by a \code{\link{MultistateModel}}
#' object, which has a \code{\link{MultistateSystem}}, which has a
#' \code{\link{TransitionMatrix}}.
#'
#' @section Fitting a model:
#' You can fit a model with 'Stan' using \code{\link{fit_stan}}.
#'
#' @section Prediction with a fitted model:
#'
#' \itemize{
#'   \item For predicting state visit risks over a time period using state path
#'  simulation, see the functions \code{\link{generate_paths}},
#'  \code{\link{p_state_visit}}, and \code{\link{p_state_visit_per_subject}}.
#'   \item For analytically solving state occupancy probabilities after some time,
#'   see \code{\link{p_state_occupancy}}.
#' }
#' @section Vignettes:
#' See the vignettes for more information.
#'
#' @name bmstate-package
#' @aliases bmstate
#' @import ggplot2
#' @import tibble
#' @import ggdist
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws_rvars
#' @importFrom dplyr mutate transmute summarize summarise
#' @importFrom dplyr inner_join semi_join left_join anti_join join_by
#' @importFrom dplyr dense_rank if_else add_count
#' @importFrom stringr str_extract str_detect str_wrap str_c
#' @importFrom survival strata Surv
#' @importFrom rlang .data
#' @importFrom tidyr expand_grid
#'
#'
"_PACKAGE"


#' Simulate example data
#'
#' @export
#' @param N number of subjects
#' @param beta_haz see the documentation of \code{\link{MultistateModel}}
#' @param beta_pk see the documentation of \code{\link{MultistateModel}}
#' @param w0 see the documentation of \code{\link{MultistateModel}}
#' @param ... arguments passed to \code{\link{create_msm}}
#' @param tm A \code{\link{TransitionMatrix}}
simulate_example_data <- function(N = 10, beta_haz = NULL,
                                  beta_pk = NULL, w0 = 1e-3, tm = NULL,
                                  ...) {
  covs <- c("sex", "age")
  pk_covs <- list(
    ka = "age",
    CL = "CrCL",
    V2 = c("weight", "sex")
  )
  if (is.null(tm)) {
    tm <- transmat_survival()
  }
  tmax <- 3 * 365.25
  mod <- create_msm(tm, covs, pk_covs, t_max = tmax, ...)
  dat <- mod$simulate_data(N, beta_haz, beta_pk, w0)
  list(model = mod, data = dat)
}

#' Example simulation setup
#'
#' @export
#' @param beta_age effect of age on death rate
#' @param beta_dose effect of dose on illness rate
#' @return a list containing a model and true betas
example_sim_setup_illnessdeath <- function(beta_dose = -0.5, beta_age = 0.5) {
  checkmate::assert_number(beta_dose)
  checkmate::assert_number(beta_age)
  covariates <- c("sex", "dose_amt", "age")
  tm <- transmat_illnessdeath()
  bh_true <- matrix(0, 2, 3)
  colnames(bh_true) <- covariates
  rownames(bh_true) <- tm$states[2:3]
  bh_true[2, 3] <- beta_age
  bh_true[1, 2] <- beta_dose
  list(
    model = create_msm(tm, covariates, categ_covs = c("sex")),
    beta_haz = bh_true
  )
}

#' Train-test split
#'
#' @export
#' @param dat A \code{\link{JointData}} object
#' @param p_test Proportion of test subjects
#' @return A list with test and train \code{\link{JointData}} objects
split_data <- function(dat, p_test = 0.25) {
  checkmate::assert_class(dat, "JointData")
  checkmate::assert_number(p_test, lower = 0, upper = 1)
  all_sub <- dat$paths$unique_subjects()
  N_sub <- length(all_sub)
  i_test <- sample.int(N_sub, round(N_sub * p_test))
  test_sub <- all_sub[i_test]
  train_sub <- setdiff(all_sub, test_sub)

  # Return
  list(
    test = dat$filter(test_sub),
    train = dat$filter(train_sub)
  )
}
