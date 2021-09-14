#' Merge Multiple Chains
#' @description This function is a method for class \code{poisreg}. Merge multiple MCMC chains into a unique chain when sampling with \code{nchains > 1} is used.
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}), with \code{nchains > 1}.
#'
#' @return The function returns an object of class \code{poisreg} with a single element \code{$sim}. 
#' The returned chains (elements of \code{sim}) are obtained by appending the simulated values of each independent chain, 
#' under the assumption that they all have reached the same stationary distribution.
#' @export
#'
#' @examples 
#' library(MASS) # load the data set
#' head(epil)
#' 
#' # Simulate multiple chains by setting nchains > 1
#' fit4 = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
#'                     iter = 1000, 
#'                     nchains = 4, thin = 2)
#' # fit4 contains 4 elements with simulation ($sim, $sim2, $sim3, $sim4)
#' 
#' mcmc_diagnostics(fit4)     
#' # the Gelman-Rubin diagnostics confirms convergence of the 4 
#' # independent chains to the same stationary distribution
#' 
#' fit4b = merge_sim(fit4) 
#' str(fit4b$sim)    
#' # fit 4b contains only one element $sim, of length 1500 
#' # (which is the result of concatenating the 4 simulations, after removing the first 25% 
#' # iterations as burn-in and keeping one iteration every two).

merge_sim = function(object)
{
  burnin = 1: (object$perc_burnin * nrow(object$sim$beta))
  
  th = seq(from = 1, to = nrow(object$sim$beta[-burnin,]), by = object$thin)
  
  object$sim$beta = object$sim$beta[-burnin,][th,]
  object$sim$r = object$sim$r[-burnin][th]
  object$sim$acceptance_rate = c(object$sim$acceptance_rate)
  object$sim$time = c(object$sim$time)
    
  for(i in 2:object$nchains)
  {
    name = paste0("sim", i)
    tmp = do.call("$", args = list(object, name))
    object$sim$beta = rbind(object$sim$beta, tmp$beta[-burnin,][th,])
    object$sim$r = c(object$sim$r, tmp$r[-burnin][th])
    object$sim$acceptance_rate = c(object$sim$acceptance_rate, tmp$acceptance_rate)
    object$sim$time = c(object$sim$time, tmp$time)
    object[[name]] = NULL
  }
  object$nchains = 1
  object$sim$acceptance_rate = mean(object$sim$acceptance_rate)
  return(object)
}
