#' @keywords internal
#' @import Rcpp
.sampling = function(formula, data, y, X,
                    iter, burnin,
                    prior, pars,
                    r_start, state,
                    thin, verbose, seed,
                    trunc_lambda)
{
  #### Importance sampler ###
  if(pars$method == "IS")
  {
    if(prior$type == "gaussian")
    {
      if(verbose == TRUE) cat(paste0(cat(paste0("Sampling ", iter, " iterations \n\n"))))
      if(!is.null(seed)) set.seed(seed)
      start = Sys.time()
      run <- .ISPolya(iter, 
                     y, X,
                     b = prior$b, B = prior$B,
                     r_start = r_start, 
                     beta_start = state,
                     max_dist = pars$max_dist,
                     max_r = pars$max_r,
                     trunc_lambda = trunc_lambda)
      end = Sys.time()
    }
    else if(prior$type == "horseshoe")
    {
      if(verbose == TRUE) cat(paste0(cat(paste0("Sampling ", iter, " iterations \n\n"))))
      start = Sys.time()
      run <- .ISPolya_horseshoe(iter, 
                               y, X,
                               b = prior$b, 
                               r_start = r_start, 
                               beta_start = state,
                               max_dist = pars$max_dist,
                               max_r = pars$max_r,
                               trunc_lambda = trunc_lambda,
                               tau = prior$tau)
      end = Sys.time()
      prior$B = NULL
    }
    
    w = exp( run$logw - mean(run$logw[-(1:round(iter/10))]) )
    run$w = w / mean(w)
    run$beta = apply(run$beta, 2, function(x) x * run$w)

    run$time = end-start
  }
  
  
  
  
  ### Metropolis-Hastings ###
  
  # gaussian
  if( (pars$method == "MH") & (prior$type == "gaussian") )
  {
    if( is.null(burnin) )
    {
      if(verbose == TRUE) cat(paste0(cat(paste0("Sampling ", iter, " iterations \n\n"))))
      if(!is.null(seed)) set.seed(seed)
      start = Sys.time()
      run <- .MHPolya_adaptr(iter, 
                            y, X,
                            b = prior$b, B = prior$B,
                            r_start = r_start, 
                            beta_start = state,
                            max_dist = pars$max_dist,
                            max_r = pars$max_r,
                            trunc_lambda = trunc_lambda)
      end = Sys.time()
    }
    else
    {
      if( !is.numeric(burnin) | (burnin < 1) | (burnin > iter) | (round(burnin)!=burnin) ) {stop("parameter burnin must be a positive integer < iter"); error = 1}
      
      if(verbose == TRUE) cat(paste0(cat(paste0("Sampling ", iter, " iterations \n\n"))))
      if(!is.null(seed)) set.seed(seed)
      start = Sys.time()
      run <- .MHPolya_adaptr_burnin(iter, 
                                    y, X,
                                    b = prior$b, B = prior$B,
                                    r_start = r_start, 
                                    beta_start = state,
                                    max_dist = pars$max_dist,
                                    burnin = burnin,
                                    max_dist_burnin = pars$max_dist_burnin,
                                    max_r = pars$max_r,
                                    trunc_lambda = trunc_lambda)
      end = Sys.time()
    }
    run$time = end-start
  }
  
  # horseshoe
  if( (pars$method == "MH") & (prior$type == "horseshoe") )
  {
    if( is.null(burnin) )
    {
      if(verbose == TRUE) cat(paste0(cat(paste0("Sampling ", iter, " iterations \n\n"))))
      if(!is.null(seed)) set.seed(seed)
      start = Sys.time()
      run <- .MHPolya_adaptr_horseshoe(iter, 
                                      y, X,
                                      b = prior$b, 
                                      tau = prior$tau,
                                      r_start = r_start, 
                                      beta_start = state,
                                      max_dist = pars$max_dist,
                                      max_r = pars$max_r,
                                      trunc_lambda = trunc_lambda)
      end = Sys.time()
    }
    else
    {
      if( !is.numeric(burnin) | (burnin < 1) | (burnin > iter) | (round(burnin)!=burnin) ) {stop("parameter burnin must be a positive integer < iter"); error = 1}
      
      if(verbose == TRUE) cat(paste0(cat(paste0("Sampling ", iter, " iterations \n\n"))))
      if(!is.null(seed)) set.seed(seed)
      start = Sys.time()
      run <- .MHPolya_adaptr_horseshoe_burnin(iter, 
                                              y, X,
                                              b = prior$b, 
                                              tau = prior$tau,
                                              r_start = r_start, 
                                              beta_start = state,
                                              max_dist = pars$max_dist,
                                              burnin = burnin,
                                              max_dist_burnin = pars$max_dist_burnin,
                                              max_r = pars$max_r,
                                              trunc_lambda = trunc_lambda)
      end = Sys.time()
    }
    prior$B = NULL
    run$time = end-start
  }
  
  return(run)
}