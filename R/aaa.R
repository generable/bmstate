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
