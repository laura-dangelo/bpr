#' Plot Trace and Distribution of Regression Parameters
#' 
#' @param x object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#' @param ... further arguments passed to or from other methods.
#'
#' @seealso \code{\link{sample_bpr}}
#' @export
plot.poisreg = function(x, ...)
{
  plot(x$sim$beta, ...)
}