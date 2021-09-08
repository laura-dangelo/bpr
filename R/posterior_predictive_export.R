#' Compute Posterior Predictive Distribution
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#' @param new_X (optional) a data frame in which to look for variables with which to predict. 
#'
#' @seealso \code{\link{posterior_predictive.poisreg}}
#' @export
posterior_predictive <- function(object, new_X) {
  UseMethod("posterior_predictive")
}