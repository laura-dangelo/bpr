#' Title
#'
#' @param object 
#' @param perc_burnin 
#' @param newdata 
#'
#' @return
#' @export
#'
#' @examples
posterior_predictive.poisreg = function(object, perc_burnin = NULL, newdata = NULL)
{
  if(!is.null(newdata)){
    if( ncol(newdata) != ncol(object$data$X) ) stop("invalid dimension for new data")}
    
  X = object$data$X
  linpred = apply(object$sim$beta, 1, function(beta) X %*% beta)
  ynew = apply(linpred, 2, function(lp) rpois(nrow(X), exp(lp)))
  
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
  
  
  meany = rpois(nrow(X), exp( X %*% colMeans(object$sim$beta[-burnin,])) )
  
  CPO = sapply(1:length(object$data$y), function(i) .CPO_i(object$data$y[i], linpred[i,]))
  LPML = sum( log(CPO) )
  
  y_pred = NULL
  y_meanpred = NULL
  if(!is.null(newdata)) 
  {
    linpred = apply(object$sim$beta, 1, function(beta) newdata %*% beta)
    y_pred = apply(linpred, 2, function(lp) rpois(nrow(newdata), exp(lp)))
    y_meanpred = rpois(nrow(X), exp( newdata %*% colMeans(object$sim$beta[-burnin,])) )
  }

  
  out = list("data" = list("X" = object$data$X, "y" = object$data$y),
             "y_pred" = ynew,
             "y_MAP_pred" = meany,
             "diagnostics" = list("CPO" = CPO, "LPML" = LPML),
             "newdata" = list( "new_X" = newdata, "y_newdata" = y_pred,
                                "y_MAP_newdata" = y_meanpred),
             "burnin" = object$burnin)
    
  class(out) = "posterior_check"
  return(out)
}



.CPO_i = function(y, linpred_i)
{
  1/ ( sum( 1/dpois(y, lambda = exp(linpred_i), log = FALSE) )/length(linpred_i) )
}