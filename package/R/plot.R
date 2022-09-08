#' Plot Trace and Distribution of Regression Parameters
#' 
#' @param x object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#' @param ... further arguments passed to or from other methods.
#'
#' @return 
#' The function calls \code{\link[coda]{plot.mcmc}} on the matrix of sampled regression coefficients, and returns the trace of the sampled outputs and a density estimate for each variable in the chain.
#' 
#' @seealso \code{\link{sample_bpr}}, \code{\link[coda]{plot.mcmc}}
#' @export
plot.poisreg = function(x, ...)
{
  plot(x$sim$beta, ...)
}