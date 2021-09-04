#' Merge Multiple Chains
#' @description This function is a method for class \code{poisreg}. Merge multiple MCMC chains into a unique chain.
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}), with \code{nchains > 1}.
#' @param perc_burnin (optional) percentage of each chain to be discarded as burn-in.
#'
#' @return The function returns an object of class \code{poisreg} with a single element \code{$sim}. 
#' The returned chains (elements of \code{sim}) are obtained by appending the simulated values of each independent chain, 
#' under the assumption that they all have reached the same stationary distribution.
#' @export
#'
#' @examples ## For examples, see example(sample_bpr).
merge.poisreg = function(object, perc_burnin = NULL)
{
  if(!is.null(object$burnin) & is.null(perc_burnin)) { 
    burnin = 1:(object$burnin+1)
  } else if(is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:round(nrow(object$sim$beta)*perc_burnin)
  } else if(is.null(object$burnin) & is.null(perc_burnin)) {
    perc_burnin = 0.25
    burnin = 1:round(nrow(object$sim$beta)*perc_burnin)
  } else if(!is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:max( round(nrow(object$sim$beta)*perc_burnin), object$burnin)
  }
  
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
