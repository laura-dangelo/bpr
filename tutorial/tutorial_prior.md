## Prior distribution
Two different prior distributions on the regression parameters are implemented: a Gaussian prior and a horseshoe prior.

The type of prior and its parameters are selected thrugh the `prior` argument: it is a named list with elements:
* `type` : whether a `"gaussian"` or a `"horseshoe"` prior should be used;
* `b`, `B` : if a Gaussian prior is used, the mean vector and covariance matrix
* `tau` : if a horseshoe prior is used, the global shrinkage parameter.

```r
## example with a Gaussian prior
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
                  iter = 1000,
                  prior = list( type = "gaussian", b = rep(0, 6), B = diag(6) * 3 ) )
                  
## example with a horseshoe prior
fit = sample_bpr( y ~  lbase*trt + lage + V4, data = epil, 
                  iter = 1000,
                  prior = list( type = "horseshoe", tau = 1 ))
```
