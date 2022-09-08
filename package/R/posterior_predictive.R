#' Compute Posterior Predictive Distribution
#' 
#' @description This function is a method for class \code{poisreg}. Compute the posterior predictive distribution and summary statistics for 
#' posterior check of the model; 
#' optionally, it also computes
#' the predictive distribution with new values of the explanatory variables.
#'
#' @param object object of class "\code{poisreg}" (usually, the result of a call to \code{\link{sample_bpr}}).
#' @param new_X (optional) a data frame in which to look for variables with which to predict. 
#'
#' @return The call to this function returns an object of S3 class \code{posterior_check}. The object is a list with the following elements:
#' @returns \code{data} : the component from \code{object} (list with covariates \code{X} and response variable \code{y}).
#' @returns \code{y_pred} : matrix of dimension \code{[n, iter]} (with \code{n} sample size), each column is a draw from the posterior predictive distribution.
#' @returns \code{y_MAP_pred} : vector of length \code{n} containing a draw from the posterior distribution obtained using the maximum a posteriori estimates (MAP) of the parameters.
#' @returns \code{diagnostics} : list containing 2 elements: \code{CPO}, i.e. the Conditional Predictive Ordinate (Gelfand et al. 1992); and \code{LPML}, i.e. 
#' the logarithm of the pseudo-marginal likelihood (Ibrahim et al. 2014).
#' @returns \code{newdata} : if the matrix \code{new_X} of new values of the covariates is provided, list of three elements: \itemize{
#' \item{\code{new_X} : the provided matrix of explanatory variables; }
#' \item{\code{y_newdata} : a matrix of dimension \code{[nrow(new_X), iter]}, each column is a draw from the posterior predictive distribution using \code{new_X};}
#' \item{\code{y_MAP_newdata} : vector of length \code{nrow(new_X)} containing a draw from the posterior distribution obtained using the MAP estimate of the parameters, 
#' computed on the new data \code{new_X}.} }
#' @returns \code{perc_burnin} : the component from \code{object}.
#' 
#' @references 
#' Gelfand, A., Dey, D. and Chang, H. (1992), Model determination using predictive distributions with implementation via sampling-based-methods (with discussion), 
#' in ‘Bayesian Statistics 4’, University Press. \cr\cr
#' Ibrahim, J. G., Chen, M.H. and Sinha, D. (2014), Bayesian Survival Analysis, American Cancer Society.
#' 
#' @export
posterior_predictive = function(object, new_X = NULL)
{
  if(!is.null(new_X)){
    if( ncol(new_X) != ncol(object$data$X) ) stop("invalid dimension for new covariate matrix")}
    
  X = object$data$X
  linpred = apply(object$sim$beta, 1, function(beta) X %*% beta)
  ynew = apply(linpred, 2, function(lp) rpois(nrow(X), exp(lp)))
  
  burnin = 1: (object$perc_burnin * nrow(object$sim$beta))
  meany = rpois(nrow(X), exp( X %*% colMeans(object$sim$beta[-burnin,])) )
  
  CPO = sapply(1:length(object$data$y), function(i) .CPO_i(object$data$y[i], linpred[i,]))
  LPML = sum( log(CPO) )
  
  y_pred = NULL
  y_meanpred = NULL
  if(!is.null(new_X)) 
  {
    linpred = apply(object$sim$beta, 1, function(beta) new_X %*% beta)
    y_pred = apply(linpred, 2, function(lp) rpois(nrow(new_X), exp(lp)))
    y_meanpred = rpois(nrow(X), exp( new_X %*% colMeans(object$sim$beta[-burnin,])) )
  }

  
  out = list("data" = list("X" = object$data$X, "y" = object$data$y),
             "y_pred" = ynew,
             "y_MAP_pred" = meany,
             "diagnostics" = list("CPO" = CPO, "LPML" = LPML),
             "newdata" = list( "new_X" = new_X, "y_newdata" = y_pred,
                                "y_MAP_newdata" = y_meanpred),
             "perc_burnin" = object$perc_burnin)
    
  class(out) = "posterior_check"
  return(out)
}

.CPO_i = function(y, linpred_i)
{
  1/ ( sum( 1/dpois(y, lambda = exp(linpred_i), log = FALSE) )/length(linpred_i) )
}