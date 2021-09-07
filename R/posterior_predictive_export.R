#' @export
posterior_predictive <- function(object, new_X) {
  UseMethod("posterior_predictive")
}