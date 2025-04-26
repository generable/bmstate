# MAIN DOCUMENTATION PAGE -------------------------------------------------

#' The 'pkmstate' package.
#'
#' @description Bayesian multistate modeling using 'Stan'.
#' @author Juho Timonen (first.last at iki.fi) and Jacki Buros
#' @keywords multistate Stan Bayesian
#'
#' @section Getting started:
#' See the following \code{R6} classes.
#' \itemize{
#'  \item \code{\link{PathData}}: Class to hold generated paths.
#' }
#'
#'
#' @name pkmstate-package
#' @aliases pkmstate
#' @import ggplot2 tibble
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws_rvars
#' @importFrom dplyr filter mutate transmute group_by select slice count
#' @importFrom dplyr summarize summarise arrange first last rename
#' @importFrom dplyr inner_join semi_join left_join join_by
#' @importFrom dplyr distinct ungroup
#' @importFrom purrr map
#' @importFrom stringr str_extract str_detect
#'
#'
"_PACKAGE"
