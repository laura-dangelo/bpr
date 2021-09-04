#' Fitting Bayesian Poisson Regression
#' @description \code{sample_bpr} is used to generate draws from the posterior distribution of the coefficients of Poisson regression models. 
#' The method allows for Gaussian and horseshoe (Carvalho et al, 2010) prior distributions.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted. 
#' @param data data frame or matrix containing the variables in the model.
#' @param y,X response variable and matrix of covariates. Alternatively to providing formula and data, one can pass y and X.
#' @param iter number of algorithm iterations.
#' @param burnin (optional) a positive integer specifying the length of the burn-in. 
#' If a value > 1 is provided, the first \code{burnin} iterations use a different tuning parameter in order to better explore the parameters space.
#' @param prior a named list of parameters to select prior type and parameters. See the details for the \code{prior} arguments.
#' @param pars a named list of parameters to select algorithm type and tuning parameters. See the details for the \code{pars} arguments.
#' @param state optional vector providing the starting points for the chains.
#' @param thin a positive integer specifying the period for saving samples. The default is 1.
#' @param verbose logical (default = TRUE) indicating whether to print messages on the progress of the algorithm and possible convergence issues.
#' @param seed (optional) positive integer: the seed of random number generator.
#' @param nchains (optional) positive integer specifying the number of Markov chains. The default is 1.
#'
#' @return An object of S3 class \code{poisreg} containing the results of the sampling. \cr\cr
#' The function \code{\link{summary}} can be used to print a summary of the posterior inference for the regression parameters and of the algorithm diagnostics. \cr
#' The function \code{\link{mcmc_diagnostics}} prints convergence diagnostics for the sampled chains.  \cr
#' The function \code{plot} prints the trace of the sampled values and a density estimate of the regression coefficients. See \code{\link[coda]{plot.mcmc}}.\cr
#' The function \code{\link{posterior_predictive}} allows to compute the posterior predictive distributions to check the model. See also the related function \code{\link{plot.ppc}}.
#' 
#' \code{poisreg} is a list containing at least the following elements:
#' @returns \code{sim} : list of the results of the sampling. It contains the following elements: \itemize{
#' \item{ \code{beta} : \code{\link[coda]{mcmc}} object of posterior draws of the regression coefficients.}
#' \item{ \code{r} : the sequence of adaptive tuning parameters used in each iteration. }
#' \item{ \code{time} : the total amount of time to perform the simulation. }
#' }
#' @returns \code{formula}  : the \code{formula} object used.
#' @returns \code{data}  : list of the data (vector of counts and covariate matrix).
#' @returns \code{state}  : the starting points of the chain.
#' @returns \code{burnin}  : length of the used \code{burnin}.
#' @returns \code{prior}  : whether a Gaussian or horseshoe prior was used.
#' @returns \code{prior_pars}  : prior parameters.
#' @returns \code{thin}  : thinning frequency passed to the \code{thin} parameter.
#' @returns \code{nchains}  : number of chains. If \code{nchains} was chosen >1, the list will also include additional numbered \code{sim} elements, one for each sampled chain.
#' 
#' 
#' @export
#'
#' @examples 
#' # Minimal call
#' out = sample_bpr(y = y, X = X, iter = 2000)
#' summary(out)   # summary of posterior inference
#' mcmc_diagnostics(out)   # summary of MCMC convergence diagnostics
#' plot(out)   # plot the traceplot and posterior density of the regression coefficients
#' 
#' @references 
#' Carvalho, C., Polson, N., & Scott, J. (2010). The horseshoe estimator for sparse signals. Biometrika, 97(2), 465-480.
#' 
#' @importFrom coda as.mcmc
#' @import Rcpp
#' @import RcppArmadillo
#' 
#' @useDynLib bpr
sample_bpr = function(formula = NULL, data = NULL, y = NULL, X = NULL,
                          iter, burnin = NULL,
                          #
                          # prior
                          prior = list(type = "gaussian",
                                       b = NULL, B = NULL,
                                       tau = NULL),
                          #
                          # pars
                          pars = list(method = "MH", 
                                      max_dist = 50,
                                      max_r = NULL, 
                                      max_dist_burnin = 1e+6),
                          #
                          # init
                          state = NULL,
                          thin = 1,
                          verbose = TRUE,
                          seed = NULL,
                          nchains = 1)
{
  if( (is.null(y) | is.null(X)) & (is.null(formula) | is.null(data)) ) {stop("please provide either y and X or a formula and data"); error = 1}
  
  ## extract column names
  if(!is.null(X))
  {
    if( !is.null(attr( X, "dimnames" )) ) { colnamesX = attr( X, "dimnames" )[[2]] }
    else { colnamesX = paste0("var", 1:ncol(X)) }
    formula = formula(y ~ X)
  }
  
  if(is.null(y) & is.null(X) & !is.null(formula) & !is.null(data))
  {
    X = as.matrix(model.matrix(formula, data))
    aa = stats::formula(formula)[[2]]
    indy = which(attr(data, "names") == aa)
    y = data[, indy]
    colnamesX = attr( model.matrix(formula, data), "dimnames" )[[2]]
  }
  
  error = 0
  trunc_lambda = 500
  n = length(y)
  p = ncol(X)
  r_start = 50
  
  ## check for call errors
  if( is.null(prior$type) ) prior$type = "gaussian"
  if( is.null(prior$b) ) prior$b = rep(0,p)
  if( is.null(prior$B) ) prior$B = diag(p)*2
  
  if( is.null(pars$method) ) pars$method = "MH"
  if( is.null(pars$max_r) ) pars$max_r = 1e+7
  if( is.null(pars$max_dist) ) pars$max_dist = 50
  if( is.null(pars$max_dist_burnin) ) pars$max_dist_burnin = 1e+6
  
  if( !is.numeric(iter) | iter < 1 | (round(iter)!=iter) ) {stop("parameter iter must be a positive integer"); error = 1}
  if( (pars$method != "MH") & (pars$method != "IS") ) {stop("method must be either 'MH' or 'IS'"); error = 1}
  if( (prior$type != "gaussian") & (prior$type != "horseshoe") ) {stop("method must be either 'gaussian' or 'horseshoe'"); error = 1}
  if( (prior$type == "horseshoe") & (is.null(prior$tau)))  {stop("please provide a value for tau"); error = 1}
  
  if(is.null(pars$max_dist)){error = 1}
  if(error==1) return(0)
  
  ## messages if verbose = T
  beta_message = " the provided initial points."
  if( is.null(state) )
  {
    if(p < (n-5)) { state = c(summary(glm(y ~ X - 1, family = "poisson"(link = "log")))$coef[,1])
    beta_message = " the maximum likelihood estimates."}
    else {state = runif(p, -.5, .5)
    beta_message = " at random initial points."}
  }
  
  if((error == 0) & (verbose == TRUE) )
  {
    cat( paste0("Running ", pars$method, " sampler with a ", prior$type, " prior distribution.", "\n", 
                "Chains initialized at", beta_message, "\n\n" ))
  }
  
  
  #---------------#       Sampling here      #---------------#
  #----------------------------------------------------------#
  sim = sampling(formula, data, y, X,
                  iter, burnin,
                  prior, pars,
                  r_start, state,
                  thin, verbose, seed,
                  trunc_lambda)
  if(verbose == TRUE) cat("Sampling completed in", round(sim$time, 5), attr(sim$time, "units"), "\n\n")
  
  sim$beta = as.mcmc(sim$beta)
  colnames(sim$beta) = colnamesX
  #class(sim) <- "poisreg"
  
  tmp = sim$acceptance_rate
  if(!is.null(burnin)) tmp2 = sim$acceptance_rate_burnin
  sim$acceptance_rate = tmp
  if(!is.null(burnin)) sim$acceptance_rate_burnin = tmp2
  
  #----------------------------------------------------------#
  
  
  run = list()
  run$sim = sim
  run$formula = formula
  run$data = list(y = y, X=X)
  run$state = state
  run$burnin = burnin
  run$prior = prior$type
  run$prior_pars = list("b" = prior$b, "B" = prior$B, "tau" = prior$tau)
  run$thin = thin
  run$nchains = nchains
  
  
  if(pars$method == "IS") 
  {
    run$method = "Importance Sampler"
    sim$effSize <- sum(sim$w)^2 / sum(sim$w^2)
    if(sim$effSize < (iter/500)) {cat("Effective sample size is low, try changing max_dist or switching to MH to improve sampling\n")}
  }
  
  if(pars$method == "MH")
  {
    run$method = "Metropolis-Hastings"
    if(sim$acceptance_rate < 0.001) {cat("Acceptance rate is low, try changing max_dist or setting different starting points\n")}
  }
  

  
  #----------------------------------------------------------#
  
  if(nchains > 1)
  {
    if(run$method == "Importance Sampler") { warning("multiple chains not implemented for IS") } else {
      
    countsim = 2
    countloop = 1
    while(!( (countsim > nchains) | (countloop > nchains + 50) ) )
    {
      name = paste0("sim", countsim)
      
      init = state + runif(p, -2, 2)
      tmp = sampling(formula, data, y, X,
                     iter, burnin,
                     prior, pars,
                     r_start, init,
                     thin, verbose = F, seed,
                     trunc_lambda)
      
      tmp$beta = as.mcmc(tmp$beta)
      colnames(tmp$beta) = colnamesX
      #class(tmp) <- "poisreg"
      
      tmpt = tmp$acceptance_rate
      if(!is.null(burnin)) tmp2 = tmp$acceptance_rate_burnin
      tmp$acceptance_rate = tmpt
      if(!is.null(burnin)) tmp$acceptance_rate_burnin = tmp2
      
      if(tmp$acceptance_rate > 0.01) 
      {
        if(verbose == TRUE) cat("Sampling chain", countsim, "completed \n")
        
        countsim = countsim + 1
        run[[name]] = tmp
      }
      countloop = countloop + 1
    }}
    if(run$nchains > countsim -1) warning("convergence not reached for some chains")
    run$nchains = countsim -1
  }
    
    
  
  class(run) = "poisreg"
  
  return(run)
  
}