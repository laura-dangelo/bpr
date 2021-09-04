#' Print Summary of the output of sample_bpr()
#'
#' @param object 
#'
#' @return
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
    text_burnin = paste0("discarding the first", max(burnin), "% iterations ")
  } else if(is.null(object$burnin) & is.null(perc_burnin)) {
    perc_burnin = 0.25
    burnin = 1:round(nrow(object$sim$beta)*perc_burnin)
    text_burnin = paste0("discarding the first", max(burnin), "% iterations ")
  } else if(!is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:max( round(nrow(object$sim$beta)*perc_burnin), object$burnin)
    text_burnin = paste0("discarding the first", max(burnin), "% iterations ")
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
    cat(paste0("(Local shrinkage parameter psi2 for horseshoe prior estimated equal to ", round(mean(object$sim$psi2[-burnin][th]),3), ")\n\n"))
  }
  
  if(object$method == "Metropolis-Hastings") 
  {
    cat("Algorithm: \n")
    if(is.null(object$burnin)) burnin = 1
    cat(paste0(" Acceptance rate in ", nrow(object$sim$beta) - max(burnin)+1, " iterations is ", round(object$sim$acceptance_rate,4), "\n") )
    if(is.null(object$burnin) | object$thin>1) cat(paste0(" Posterior estimates computed ", text_burnin, text_thin, "\n"))
    cat(paste0(" Mean effective sample size is equal to ", round(mean(effectiveSize(object$sim$beta[-burnin,][th,])))) )
  }
  
  if(object$method == "Importance Sampler") 
  {
    cat("Algorithm: \n")
    if(is.null(object$burnin)) burnin = 1
    cat(paste0(" Posterior estimates computed discarding " , text_burnin, text_thin, "\n"))
    cat(paste0(" Mean effective sample size is equal to ", round(sum(object$sim$w[-burnin][th])^2 / sum(object$sim$w[-burnin][th]^2) )) )
  }
  
  
}


