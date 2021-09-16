# MCMC parameters
Most MCMC tuning quantities are passed through the `pars` paramter: it is a named list with elements
* `method` : used to choose the algorithm (default is Metropolis-Hastings). See [tutorial_IS](/tutorial_IS.md) for the importance sampler;
* `max_dist` : main tuning parameter, sets the "distance" of the approximation from the true target posterior;
* `max_r` : additional tuning parameter, used to additionally "further" the approximation when max_dist is not sufficient;
* `max_dist_burnin` : tuning parameter used during the burn-in to explore the parameter space. Setting a large distance favors a better exploration of the space, while a smaller distance allows to obtain a better mixing.

Other optional parameters, not passed through the list, are
* `burnin` : length of the burnin. The first part of the chain is sampled using the `max_dist_burnin` parameter.
* `thin` : frequency of the thinning.
* `perc_burnin` : percentage of the chain to discard to compute the estimates, can be different from `burnin`.

### Tuning parameter
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
