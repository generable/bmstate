#' Creates the Stan model
#'
#' @param mixp use mixture prior?
#' @export
#' @return A CmdstanR model
create_stan_model <- function(mixp = FALSE) {
  fn <- "msm-general.stan"
  if (mixp) {
    fn <- "msm-general-mixp.stan"
  }
  filepath <- system.file(file.path("stan", fn), package = "pkmstate")
  cmdstanr::cmdstan_model(filepath)
}
