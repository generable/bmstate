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
#'  \item \code{\link{PathData}}
#'  \item \code{\link{create_brier_score_plot}}
#'  \item \code{\link{create_cindex_plot}}
#'  \item \code{\link{create_stan_data}}
#'  \item \code{\link{create_stan_model}}
#'  \item \code{\link{do_split}}
#'  \item \code{\link{fit_coxph}}
#'  \item \code{\link{fit_mstate}}
#'  \item \code{\link{generate_paths_many_subjects}}
#'  \item \code{\link{legend_to_PT_matrix}}
#'  \item \code{\link{legend_to_TFI_matrix}}
#'  \item \code{\link{plot_cumhaz_msfit}}
#'  \item \code{\link{plot_other_beta}}
#'  \item \code{\link{plot_p_event}}
#'  \item \code{\link{plot_p_event_by_dose}}
#'  \item \code{\link{simulate_example_data}}
#'  \item \code{\link{subject_df_with_idx}}
#'  \item \code{\link{summarize_ppsurv}}
#' }
#'
#' @name bmstate-package
#' @aliases bmstate
#' @import ggplot2
#' @import tibble
#' @import ggdist
#' @importFrom ggpubr ggarrange get_legend annotate_figure
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws_rvars
#' @importFrom dplyr filter mutate transmute group_by select slice count
#' @importFrom dplyr summarize summarise arrange first last rename
#' @importFrom dplyr inner_join semi_join left_join anti_join join_by
#' @importFrom dplyr distinct ungroup dense_rank if_else bind_rows
#' @importFrom dplyr count add_count pull lead
#' @importFrom purrr map
#' @importFrom stringr str_extract str_detect str_wrap str_c
#' @importFrom survival strata
#' @importFrom rlang .data
#' @importFrom tidyr expand_grid
#'
#'
"_PACKAGE"
