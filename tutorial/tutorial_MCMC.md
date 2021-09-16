# MCMC parameters
Most MCMC tuning quantities are passed through the `pars` paramter: it is a named list with elements
* `method` : used to choose the algorithm (default is Metropolis-Hastings). See [tutorial_IS](/tutorial_IS.md) for the importance sampler;
* `max_dist` : main tuning parameter, sets the "distance" of the approximation from the true target posterior;
* `max_r` : additional tuning parameter, used to additionally worsen the approximation when max_dist is not sufficient;
* `max_dist_burnin` : tuning parameter used during the burn-in to explore the parameter space. Setting a large distance favors a better exploration of the space, while a smaller distance allows to obtain a better mixing.

Other optional parameters, not passed through that list, are
* `burnin` : length of the burn-in. The first part of the chain is sampled using the `max_dist_burnin` parameter;
* `thin` : frequency of the thinning;
* `perc_burnin` : percentage of the chain to discard to compute the estimates, can be different from `burnin`;
* `nchains` : number of MCMC chains.

## Tuning parameters
The main tuning parameter of the algorithm is `max_dist`, which allows to balance acceptance rate and autocorrelation to obtain a good mixing of the chains.
```r
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, iter = 1000,
                   pars = list( max_dist = 10 ))
fit$sim$acceptance_rate
## [1] 0.489
```

To have a higher acceptance rate, `max_dist` must be small: in this case, the algorithm is an approximate Gibbs sampler. However, autocorrelation is usually higher.
```r
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, iter = 1000,
                   pars = list( max_dist = 0.01 ))
fit$sim$acceptance_rate
## [1] 0.997
```

To have, on the contrary, a low acceptance rate, sometimes only setting a large distance is not enough, in this case, the parameter `max_r` allows to keep the approximation very far from the target posterior:
```r
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, iter = 1000,
                   pars = list( max_dist = 5000 ))
fit$sim$acceptance_rate
## [1] 0.443
   
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, iter = 1000,
                   pars = list( max_dist = 5000, max_r = 10 ))
fit$sim$acceptance_rate
## [1] 0.059
```

## Burn-in
When the `burnin` option is used, a different tuning parameter is used for the first part (the burn-in) and the last part (the actual sampling) of the chain. In particular, for the first part it is used the parameter `max_dist_burnin`, while for the second part is used `max_dist`.

A different parameter is `perc_burnin`, which only defines the proportion of the chain to discard when computing the estimates.

```r
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
                  iter = 1000,
                  burnin = 200, # use max_dist_burnin for the first 200 iterations
                  pars = list( max_dist = 1, max_dist_burnin = 1e+7 ),
                  perc_burnin = 0.5 # discard the first 50% of the chain to do inference
                  )
```


## Multiple chains
To run multiple chains and see whether they all converged to the same stationary distribution, it is possible to use the `nchains` parameter:
```r
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
                  iter = 1000,
                  pars = list( max_dist = 10 ),
                  nchains = 4)
```
In this case, the MCMC summary also prints an additional section with the Gelman-Rubin test of convergence:
```r
mcmc_diagnostics(fit)

   ## Total chains length = 1000 
   ## Discarding the first 250 iterations as burnin 
   ## Thinning frequency = 1 
   ## 
   ## MCMC Diagnostics: 
   ##                    Eff. Size Geweke test Pr(>|z|)
   ## (Intercept)           214.13      -0.189   0.8500
   ## lbase                 200.83       0.101   0.9196
   ## trtprogabide          172.88       1.218   0.2234
   ## lage                  119.67      -1.062   0.2884
   ## V4                    193.04       2.349   0.0188
   ## lbase:trtprogabide    179.43      -3.119   0.0018
   ## 
   ## 
   ## Gelman-Rubin convergence diagnostic on 4 chains 
   ##  Potential scale reduction factors:
   ## 
   ##                    Point est. Upper C.I.
   ## (Intercept)              1.00       1.01
   ## lbase                    1.00       1.01
   ## trtprogabide             1.00       1.01
   ## lage                     1.00       1.00
   ## V4                       1.01       1.02
   ## lbase:trtprogabide       1.00       1.00
   ## 
   ## Multivariate psrf
   ## 
   ## 1.01
```

Once the convergence to a common distribution has been assessed, it is possible to automatically merge the 4 chains by calling `merge_sim(fit)`.
















