## `bpr` tutorial

This tutorial is a short guide on the use of the `bpr` package to perform posterior sampling and inference on the parameters of Bayesian Poisson regression models.

### Model specification
First, we briefly recall the model specification and prior distributions available in the package.
The general setting is that of standard Poisson regression models, where the interest is to regress a vector of counts *y*, of length *n*, on an *n x p* matrix *X* of covariates:

<p align="center">
y ~ Poisson(β)

λ = exp{ X * β }
</p>

where β is a length-*p* vector of regression coefficients, and λ is the linear predictor (vector of length *n*).

On the regression coefficients we place conditionally Gaussian prior distributions: we assume that, conditionally on the mean vector *b* and covariance matrix *B*,
<p align="center">
(β | b, B) ~ N(b,B).
</p>

Implemented in this package are a straightforward Gaussian prior distribution with informative *(b,B)* fixed using prior information, and the horseshoe prior of Carvalho et al. (2010).
The horseshoe prior is a scale mixture of Gaussians where *b* is set to zero and the variance has a hierarchical representation:
<p align="center">
(β<sub>j</sub> | η<sup>2</sup><sub>j</sub>,τ<sup>2</sup>) ~ N(0,η<sup>2</sup>τ<sup>2</sup>)

η ~ C+(0,1),  τ ~ C+(0,1)
</p>





#### References
Carvalho, C. M., Polson, N. G. and Scott, J. G. (2010), "The horseshoe estimator for sparse signals", *Biometrika* **97**(2), 465–480
