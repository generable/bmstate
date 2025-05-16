# MAIN DOCUMENTATION PAGE -------------------------------------------------
#
# dev note: see
#   https://cran.r-project.org/web/packages/dplyr/vignettes/in-packages.html
# for how to use dplyr in a package

#' The 'bmstate' package.
#'
#' @description Bayesian multistate modeling using 'Stan'.
#' @author Juho Timonen (first.last at iki.fi) and Jacki Buros
#' @keywords multistate Stan Bayesian
#'
#' @section Getting started:
#' See the vignettes. You can create a model using \code{\link{create_msm}}.
#' The following classes and their documentation can be useful
#'
#' \itemize{
#'  \item \code{\link{MultistateModel}}
#'  \item \code{\link{PKModel}}
#'  \item \code{\link{PathData}}
#'  \item \code{\link{TransitionMatrix}}.
#' }
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
simulate_example_data <- function(N = 10, beta_haz = NULL,
                                  beta_pk = NULL, w0 = 1e-3) {
  covs <- c("sex", "age")
  pk_covs <- list(
    ka = "age",
    CL = "CrCL",
    V2 = c("weight", "sex")
  )
  tm <- transmat_full(self_loops = FALSE)
  tmax <- 3 * 365.25
  mod <- create_msm(tm, covs, pk_covs, FALSE, tmax = tmax)
  mod$simulate_data(N, beta_haz, beta_pk, w0)
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
    model = create_msm(tm, covariates),
    beta_haz = bh_true
  )
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

#' Plot event time distribution
#'
#' @export
#' @param t a vector of event times
plot_time_dist <- function(t) {
  checkmate::assert_numeric(t)
  ggplot(data.frame(Time = t), aes(x = .data$Time)) +
    ggdist::stat_halfeye() +
    geom_vline(
      mapping = NULL,
      xintercept = mean(t),
      color = "gray20",
      lty = 2
    )
}


#' Create the main Stan model
#'
#' @export
#' @param ... Arguments passed to CmdStanR
create_stanmodel <- function(...) {
  filepath <- system.file(file.path("stan", fn), package = "bmstate")
  # silence compile warnings from cmdstan
  utils::capture.output(
    {
      mod <- cmdstanr::cmdstan_model(filepath, ...)
    },
    type = "message"
  )
  mod
}
