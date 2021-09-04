#' MCMC Convergence Diagnostics
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#' @param perc_burnin (optional) percentage of each chain to be discarded as burn-in.
#'
#' @return \code{mcmc_diagnostics.poisreg} returns a list with elements:
#' @returns \code{chain_length} total length of the MCMC chains.
#' @returns \code{len_burnin}  the length of the burn-in used to compute the estimates.
#' @returns \code{thin} the thinning frequency used (from \code{object}).
#' @returns \code{effSize} effective sample size of each parameter chain after removing burn-in and thinning.
#' @returns \code{geweke} Geweke diagnostics of convergence of the chains (value of the test and p-value).
#' @returns \code{gelman_rubin} if \code{nchains > 1}, Gelman-Rubin diagnostics of convergence.

#' @export
#'
#' @export
#' @importFrom coda effectiveSize geweke.diag gelman.diag
mcmc_diagnostics.poisreg = function(object, perc_burnin = NULL)
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

