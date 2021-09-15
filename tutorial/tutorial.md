## `bpr` tutorial

This tutorial is a short guide on the use of the `bpr` package to perform posterior sampling and inference on the parameters of Bayesian Poisson regression models.

### Model specification
First, we briefly recall the model specification and prior distributions available in the package.
The general setting is that of standard Poisson regression models, where the interest is to regress a vector of counts *y*, of length *n*, on an *n x p* matrix *X* of covariates:

<p align="center">
y ~ Poisson( λ )
</p> <p align="center">
λ = exp{ X * β }
</p>

where β is a length-*p* vector of regression coefficients, and λ is the linear predictor (vector of length *n*).

On the regression coefficients we place conditionally Gaussian prior distributions: we assume that, conditionally on the mean vector *b* and covariance matrix *B*,
<p align="center">
( β | b, B ) ~ N( b, B ).
</p>

Implemented in this package are a straightforward Gaussian prior distribution with informative *(b,B)* fixed using prior information, and the horseshoe prior of Carvalho et al. (2010).
The horseshoe prior is a scale mixture of Gaussians where *b* is set to zero and the variance has a hierarchical representation:
<p align="center">
( β<sub>j</sub> | η<sub>j</sub><sup>2</sup>, τ<sup>2</sup> ) ~ N( 0, η<sup>2</sup> τ<sup>2</sup> )
</p><p align="center">
η ~ C<sup>+</sup>( 0, 1 )  ,                 τ ~ C<sup>+</sup>( 0, 1 )
</p>


### Example
We use the `epil` data set from the `MASS` library, containing seizure counts for 59 epileptics.
```{r}
library(bpr)
library(MASS)

head(epil)
```
   ##   y     trt base age V4 subject period      lbase       lage
   ##   5 placebo   11  31  0       1      1 -0.7563538 0.11420370
   ##   3 placebo   11  31  0       1      2 -0.7563538 0.11420370
   ##   3 placebo   11  31  0       1      3 -0.7563538 0.11420370
   ##   3 placebo   11  31  1       1      4 -0.7563538 0.11420370
   ##   3 placebo   11  30  0       2      1 -0.7563538 0.08141387
   ##   5 placebo   11  30  0       2      2 -0.7563538 0.0814138

#### Model 1: Metropolis-Hastings algorithm with Gaussian priors
A minimal call only requires to specify the data and number of MCMC iterations.
Calling the main function `sample_bpr()` with the default parameters will implement a Metropolis-Hastings algorithm with tuning parameter `max_dist` equal to 50, and with independent *N(0,2)* prior distributions on the regression parameters.
```r
fit = sample_bpr( y ~ lbase * trt + lage + V4, data = epil, 
                    iter = 1000)
```
if `verbose = TRUE` (default) the call will produce a very synthetic output:
   ## Running MH sampler with a gaussian prior distribution.
   ## Chains initialized at the maximum likelihood estimates.
   ##
   ## Sampling 1000 iterations 
   ##
   ## Sampling completed in 0.114 secs
   
If there are important issues with the algorithm performance (such as acceptance rate equal to 0), a warning message is also printed.





















### References
Carvalho, C. M., Polson, N. G. and Scott, J. G. (2010), "The horseshoe estimator for sparse signals", *Biometrika* **97**(2), 465–480
