#' Graphical Posterior Predictive Checks
#' @description This function is a method for class \code{posterior_check}. Plot diagnostic statistics for graphical posterior predictive checks.
#'
#' @param x object of class "\code{posterior_check}" (usually, the result of a call to \code{\link{posterior_predictive}}).
#' @param ... other parameters to be passed through to plotting functions. See Details.
#' 
#' 
#' @return 
#' The function outputs (at least) three plots for graphical posterior predictive check.\cr
#' The first plot compares the empirical cumulative distribution function (ECDF) with the cumulative distribution function obtained 
#' with samples from the posterior predictive distribution (median and point-wise 95\% credible bands). \cr
#' The second plot compares the distribution of the observed sample with the predictive distribution obtained using the maximum a posteriori (MAP) 
#' estimates of the regression parameters.\cr
#' The third plot compares the predictive distribution of a statistic (default is the mean) with the observed value of the same statistics,
#'  displayed with a red line.
#' 
#' 
#' @details 
#' It is possible to generate additional plots that compare the posterior predictive distribution of a statistics with the observed value. 
#' This is done through the parameter \code{stats}: it is a list with elements the function names of the statistics one wants to compare. 
#' Default is \code{stats = list("mean")}, other possible values are, e.g., "median", "sd", "max" etc.
#' 
#' @seealso 
#' \code{\link{posterior_predictive}}
#' 
#' @examples 
#' library(MASS) # load the data set
#' head(epil)
#' 
#' fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
#'                    iter = 1000)
#' \donttest{ plot(posterior_predictive(fit), stats = c("mean", "sd", "max"))   }
#' # plots for posterior predictive check
#' 
#' 
#' @export
#' @importFrom graphics abline hist legend lines mtext par polygon
#' @importFrom stats dpois ecdf glm median model.matrix pnorm poisson rpois runif sd stepfun
#' @importFrom grDevices devAskNewPage
plot.posterior_check = function(x, ...)
{
  burnin = 1: (x$perc_burnin * ncol(x$y_pred))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  .plotECDF(x, burnin)
  devAskNewPage(ask = TRUE)
  .plotHIST(x, burnin)
  .plotSTAT(x, burnin, ...)
  devAskNewPage(ask = FALSE)
}

.plotECDF = function(object, burnin)
{
  par(mfrow = c(1,1))
  seqq = sort(unique(c(0,object$y_pred[,-burnin],object$data$y)))
  ecdfseq = apply( object$y_pred[,-burnin], 2, function(yy) sapply(seqq, function(x) ecdf(yy)(x) ) ) 
  
  median_cdf = apply(ecdfseq, 1, median)
  low_cdf = apply(ecdfseq, 1, function(x) .emp.hpd(x)[1])
  upp_cdf = apply(ecdfseq, 1, function(x) .emp.hpd(x)[2])
  
  plot(ecdf(object$data$y), verticals = TRUE, col.hor = "salmon",col.ver = "salmon",
       main = "CDF\n", do.points = F,
       xlab = "y", ylab = "Fn(y)")
  subtit =  paste0("ECDF vs. posterior predictive cdf (median and pointwise 95% credible bands)")
  mtext(side = 3, line = 0.5,  subtit)
  
  for(i in 1:length(seqq))
  {
    polygon(seqq[c(i,i,i+1,i+1)], c(low_cdf[i], upp_cdf[i], upp_cdf[i], low_cdf[i]), col = "#b6d1db", border = NA)
  }

  lines(stepfun(seqq+0.07, c(0,median_cdf)), col = "darkslategrey", lwd = 2, cex = 0.7)
  lines(ecdf(object$data$y), verticals = TRUE, bg = "salmon", lwd = 2, col.points = "salmon", pch=1, cex = 0.7,
        col.hor = "salmon",col.ver = "salmon")
  
  legend("bottomright", c("observed", "predicted"), border = NA,
         lty = 1,
         lwd = 2,
         pch = 21, seg.len = 0.5,
         bg = c("salmon", "#b6d1db"),
         fill = c(NA, "#b6d1db"),
         col = c("salmon", "darkslategrey"),
         ncol = 1, bty = "n")
  
}

.plotHIST = function(object, burnin)
{
  par(mfrow = c(1,2))
  
  hist(object$data$y, main = "Observed y", xlab = "y")
  hist(object$y_MAP_pred, main = "MAP predictive distribution", xlab = "y")
}

.plotSTAT = function(object, burnin, stats = c("mean"))
{
  par(mfrow = c(1,1))
  if(length(stats) == 1){
    res = apply(object$y_pred[,-burnin], 2, function(y) do.call(stats, list(y)) )
    hist(res, main = paste0("Histogram of posterior predictive ", stats),
         xlab = stats, xlim = range(c(res, do.call(stats, list(object$data$y)))) )
    subtit =  paste0("Observed value vs. posterior predictive distribution")
    mtext(side = 3, line = 0.5,  subtit)
    abline(v = do.call(stats, list(object$data$y)), col = 2, lwd = 1.5)
    legend("topright", c("observed"), lty=1, col = 2, bty = "n")
  }
  else{
    for(i in 1:length(stats)) {
      res = apply(object$y_pred[,-burnin], 2, function(y) do.call(stats[i], list(y)) )
      hist(res, main = paste0("Histogram of posterior predictive ", stats[i] ),
           xlab = stats[i], xlim = range(c(res, do.call(stats[i], list(object$data$y)))) )
      subtit =  paste0("Observed vs. posterior predictive distribution")
      mtext(side = 3, line = 0.5,  subtit)
      abline(v = do.call(stats[i], list(object$data$y)), col = 2, lwd = 1.5)
      legend("topright", c("observed"), lty=1, col = 2, bty = "n")
    }
  }
}



 