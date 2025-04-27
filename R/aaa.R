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
#' See the vignettes and the following exported functions and classes
#'
#' \itemize{
#'  - \code{\link{PathData}}
#'  - \code{\link{create_stan_data}}
#'  - \code{\link{create_stan_model}}
#'  - \code{\link{do_split}}
#'  - \code{\link{fit_coxph}}
#'  - \code{\link{fit_mstate}}
#'  - \code{\link{generate_paths_many_subjects}}
#'  - \code{\link{legend_to_PT_matrix}}
#'  - \code{\link{legend_to_TFI_matrix}}
#'  - \code{\link{plot_cumhaz_msfit}}
#'  - \code{\link{plot_p_event}}
#'  - \code{\link{plot_p_event_by_dose}}
#'  - \code{\link{simulate_example_data}}
#'  - \code{\link{subject_df_with_idx}}
#'  - \code{\link{summarize_ppsurv}}
#' }
#'
#' @name bmstate-package
#' @aliases bmstate
#' @import ggplot2 tibble
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws_rvars
#' @importFrom dplyr filter mutate transmute group_by select slice count
#' @importFrom dplyr summarize summarise arrange first last rename
#' @importFrom dplyr inner_join semi_join left_join join_by
#' @importFrom dplyr distinct ungroup dense_rank
#' @importFrom purrr map
#' @importFrom stringr str_extract str_detect
#' @importFrom survival strata
#' @importFrom rlang .data
#'
#'
"_PACKAGE"
