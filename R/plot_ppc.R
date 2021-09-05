#' Plot posterior check
#'
#' @param object 
#' @param perc_burnin 
#' @param ... 
#' 
#' @export
plot.posterior_check = function(object, perc_burnin = NULL, ...)
{
  if(!is.null(object$burnin) & is.null(perc_burnin)) { 
    burnin = 1:(object$burnin+1)
  } else if(is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:round(ncol(object$y_pred)*perc_burnin)
  } else if(is.null(object$burnin) & is.null(perc_burnin)) {
    perc_burnin = 0.25
    burnin = 1:round(ncol(object$y_pred)*perc_burnin)
  } else if(!is.null(object$burnin) & !is.null(perc_burnin)) {
    burnin = 1:max( round(ncol(object$y_pred)*perc_burnin), object$burnin)
  }
  
  .plotECDF(object, burnin)
  devAskNewPage(ask = TRUE)
  .plotHIST(object, burnin)
  .plotSTAT(object, burnin, ...)
  devAskNewPage(ask = FALSE)
}

.plotECDF = function(object, burnin)
{
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
  par(mfrow = c(1,1))
}

.plotSTAT = function(object,  burnin, stats = c("mean"))
{
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



 