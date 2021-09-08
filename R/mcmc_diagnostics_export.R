#' MCMC Convergence Diagnostics
#' @param object an object for which mcmc diagnostic is desired.
#' @seealso \code{\link{mcmc_diagnostics.poisreg}}
#'
#' @export
mcmc_diagnostics <- function(object) {
  UseMethod("mcmc_diagnostics")
}