#' MCMC Convergence Diagnostics
#' @description This function is a method for class \code{poisreg}. It prints convergence diagnostics and accuracy statistics of the MCMC output.
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#'
#' @details The printed output of \code{mcmc_diagnostics} summarizes some common convergence diagnostics for Markov chains. 
#' The first part recaps the total length, burn-in and thinning used for the simulation.
#' 
#' The second part is a table with diagnostic statistics about each chain of the regression parameters. The first column is 
#' the effective sample size computed after removing the burn-in and thinning. 
#' The last two columns report the value and observed p-value of the Geweke test of equality of the first and last part of the chain.
#' 
#' The last part is printed only if multiple chains are computed. In this case, it reports the Gelman-Rubin statistics to test convergence to the same stationary  
#' distribution. Values much larger than 1 suggest lack of convergence to a common distribution.
#'
#' @return \code{mcmc_diagnostics} returns a list with elements:
#' @returns \code{chain_length} : total length of the MCMC chains.
#' @returns \code{len_burnin} : the length of the burn-in used to compute the estimates.
#' @returns \code{thin} : the thinning frequency used (from \code{object}).
#' @returns \code{effSize} : effective sample size of each parameter chain after removing burn-in and thinning. See \code{\link[coda]{effectiveSize}}.
#' @returns \code{geweke} : Geweke diagnostics of convergence of the chains (value of the test and p-value). See \code{\link[coda]{geweke.diag}}
#' @returns \code{gelman_rubin} : if \code{nchains > 1}, Gelman-Rubin diagnostics of convergence. See \code{\link[coda]{gelman.diag}}.
#'
#' @seealso \code{\link{summary.poisreg}} , \code{\link{plot.poisreg}} ,
#' \code{\link{merge_sim}} , \code{\link[coda]{effectiveSize}} , \code{\link[coda]{geweke.diag}} , \code{\link[coda]{gelman.diag}}
#' @examples 
#' # For examples see example(sample_bpr)
#' 
#' @export
#' @importFrom coda effectiveSize geweke.diag gelman.diag
mcmc_diagnostics = function(object)
{
  burnin = 1: (object$perc_burnin * nrow(object$sim$beta))
  th = seq(from = 1, to = nrow(object$sim$beta[-burnin,]), by = object$thin)

  if(object$method == "Metropolis-Hastings")
  {
    cat("\n")
    cat("Total chains length =", nrow(object$sim$beta), "\n")
    cat("Discarding the first", max(burnin), "iterations as burnin \n")
    cat("Thinning frequency =", object$thin, "\n\n")
    tmp <- do.call(data.frame, 
                   list( effSize = round(effectiveSize(object$sim$beta[-burnin,][th,]), 2),
                         geweke = format(c(geweke.diag(object$sim$beta[-burnin,][th,])$z), digits = 3),
                         pr = format(pnorm(abs(geweke.diag(object$sim$beta[-burnin,][th,])$z), lower.tail = F)*2, digits = 2)
                         )
    )
    colnames(tmp) = c("Eff. Size", "Geweke test", "Pr(>|z|)")
    rownames(tmp) = colnames(object$sim$beta)
    cat("MCMC Diagnostics: \n")
    print(tmp)
    cat("\n\n")
  }
  else if(object$method == "Importance Sampler")
  {
    effS = round( sum(object$w[-burnin][th])^2 / sum(object$w[-burnin][th]^2), 2 )
    cat("Effective sample size computed on the weights is equal to", effS)
  }
  
  if(object$nchains > 1)
  {
    multch = list()
    multch[[1]] = object$sim$beta
    for(i in 2:object$nchains)
    {
      name = paste0("sim", i)
      multch[[i]] = do.call("$", args = list(object, name))
      multch[[i]] = multch[[i]]$beta
    }
    cat(paste0("Gelman-Rubin convergence diagnostic on ", object$nchains, " chains \n "))
    print(gelman.diag(multch))
  }
  
  out = list()
  out$chain_length = nrow(object$sim$beta)
  out$len_burnin = max(burnin)
  out$thin = object$thin
  
  out$effSize = effectiveSize(object$sim$beta[-burnin,][th,])
  out$geweke = data.frame( "Geweke test" = c(geweke.diag(object$sim$beta[-burnin,][th,])$z),
        "Pr(>|z|)" = pnorm(abs(geweke.diag(object$sim$beta[-burnin,][th,])$z), lower.tail = F)*2  )
  
  if(object$nchains > 1) out$gelman_rubin = gelman.diag(multch)
  
  return(invisible(out))
}
