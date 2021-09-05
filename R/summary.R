#' Summarizing Bayesian Poisson Regression Fit
#' @description This function is a method for class \code{poisreg}. It prints summary statistics and returns posterior estimates of regression quantities.
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#' @param perc_burnin (optional) percentage of each chain to be discarded as burn-in. Default is 0.25.
#' 
#' @details The printed output of \code{summary.poisreg} summarizes the main quantities of the fit. 
#' The first component \code{Call} recaps the type of prior and algorithm used.
#' 
#' \code{Coefficients} is a table of estimated quantities for the regression parameters. The first three columns report the estimated posterior mean, 
#' standard errors and medians. The last two columns correspond to the lower and upper bounds of the 0.95 credible intervals. 
#' If the credible interval does not include zero, a star is printed in correspondence of each parameter 
#' (similarly to the 'significance stars' of \code{\link[stats]{summary.lm}}).
#' All the estimates are computed discarding the first part of the chain as burn-in (more details are printed in the \code{Algorithm} section).
#' 
#' \code{Algorithm} briefly summarizes the main diagnostics of convergence and efficiency of the algorithm. 
#' It prints the number of iterations actually used to obtain the estimates, after removing the burn-in and thinning.
#' If a Metropolis-Hastings algorithm is used, the summary reports the acceptance rate, 
#' which is the most commonly used indicator to tune the performance of the algorithm, along with the mean effective sample size 
#' (averaged over all parameters).
#' If the importance sampler is used, the summary only reports the effective sample size, which is a computed as \eqn{\sum_{t} w_t^2 / (\sum_{t} w_t)^2}
#' (where \eqn{w_t} is the sequence of weights) and is a measure of the efficiency of the sampler.
#' 
#' 
#' @return \code{summary.poisreg} returns a list with elements:
#' @returns \code{formula} the component from \code{object}.
#' @returns \code{data} the component from \code{object} (list with covariates \code{X} and response variable \code{y}).
#' @returns \code{prior} \code{prior$type} from \code{object}.
#' @returns \code{prior_pars} prior parameters from \code{object}.
#' @returns \code{coefficients} the matrix of coefficients, standard errors and 95\% credible intervals.
#' @returns \code{psi2} if a horseshoe prior is selected, the estimate of the local shrinkage parameter.
#' @returns \code{len_burnin} the length of the burn-in used to compute the estimates.
#' @returns \code{effSize} the mean effective sample size of the chains used to compute the estimates.
#' 
#' @examples 
#' # For examples see example(sample_bpr)
#'
#' 
#' @export
#'
#' @importFrom coda effectiveSize
summary.poisreg <- function(object, perc_burnin = NULL) {
  cat("\n")
  cat("Call: \n", "sample_poisreg( formula = ",  deparse(object$formula), 
      ", prior = ", object$prior,
      ", algorithm = ", object$method, ")", "\n\n")
  
  if(!is.null(object$burnin) & is.null(perc_burnin)) { 
    burnin = 1:(object$burnin+1)
    text_burnin = "" 
  } else if(is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:round(nrow(object$sim$beta)*perc_burnin)
    text_burnin = paste0("discarding the first ", max(burnin), " iterations ")
  } else if(is.null(object$burnin) & is.null(perc_burnin)) {
    perc_burnin = 0.25
    burnin = 1:round(nrow(object$sim$beta)*perc_burnin)
    text_burnin = paste0("discarding the first ", max(burnin), " iterations ")
  } else if(!is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:max( round(nrow(object$sim$beta)*perc_burnin), object$burnin)
    text_burnin = paste0("discarding the first ", max(burnin), " iterations ")
  }
  
  
  th = seq(from = 1, to = nrow(object$sim$beta[-burnin,]), by = object$thin)
  text_thin = ""
  if(object$thin > 1) text_thin = paste0("with thinning frequency = ", object$thin)
  
  tmp <- do.call(data.frame, 
                 list( mean = format(colMeans(object$sim$beta[-burnin,][th,]), digits = 5),
                       stderror = format(apply(object$sim$beta[-burnin,][th,], 2, sd), digits = 5),
                       median = format(apply(object$sim$beta[-burnin,][th,], 2, median), digits = 5),
                       lowCI = format(apply(object$sim$beta[-burnin,][th,], 2, .emp.hpd)[1,], digits = 3),
                       uppCI = format(apply(object$sim$beta[-burnin,][th,], 2, .emp.hpd)[2,], digits = 3),
                       iszero = rep("*", ncol(object$sim$beta)))
  )
  tmp$iszero[as.numeric(tmp$lowCI) * as.numeric(tmp$uppCI) < 0] = ""
  colnames(tmp) = c("Mean", "Std. Error", "Median", "Lower CI", "Upper CI", "")
  rownames(tmp) = colnames(object$sim$beta)

  cat("Coefficients: \n")
  print(tmp)
  cat("--- \n '*' if 95% credible interval does not include zero. \n \n")
  
  if(object$prior == "horseshoe")
  {
    psi2est = round(mean(object$sim$psi2[-burnin][th]),3)
    cat(paste0("(Local shrinkage parameter psi2 for horseshoe prior estimated equal to ", psi2est, ")\n\n"))
  }
  
  effSize = NULL
  if(object$method == "Metropolis-Hastings") 
  {
    effSize = round(mean(effectiveSize(object$sim$beta[-burnin,][th,])))
    cat("Algorithm: \n")
    if(is.null(object$burnin)) burnin = 1
    if(is.null(object$burnin) | object$thin>1) cat(paste0(" Posterior estimates computed ", text_burnin, text_thin, "\n"))
    cat(paste0(" Acceptance rate in ", nrow(object$sim$beta) - max(burnin)+1, " iterations is ", round(object$sim$acceptance_rate,4), "\n") )
    cat(paste0(" Mean effective sample size is equal to ", effSize) )
  }
  
  if(object$method == "Importance Sampler") 
  {
    effSize = round(sum(object$sim$w[-burnin][th])^2 / sum(object$sim$w[-burnin][th]^2) )
    cat("Algorithm: \n")
    if(is.null(object$burnin)) burnin = 1
    cat(paste0(" Posterior estimates computed discarding " , text_burnin, text_thin, "\n"))
    cat(paste0(" Mean effective sample size is equal to ", effSize) )
  }
  
  out = list()
  out$formula = object$formula
  out$data = object$data
  out$prior = object$prior
  out$prior_pars = object$prior_pars
  
  out$coefficients = tmp
  if(object$prior == "horseshoe") out$psi2 = psi2est
  out$len_burnin = max(burnin)
  out$effSize = effSize
  
  return(invisible(out))
}


