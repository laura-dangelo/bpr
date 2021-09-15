## `bpr` tutorial

This tutorial is a short guide on the use of the `bpr` package to perform posterior sampling and inference on the parameters of Bayesian Poisson regression models.

### Model specification
First, we briefly recall the model specification and prior distributions available in the package.
The general setting is that of standard Poisson regression models, where the interest is to regress a vector of counts *y*, of length *n*, on an *n x p* matrix *X* of covariates:

<p align="center">
y ~ Poisson(β)
lambda = exp{ X * β }
</p>

where β is a length-*p* vector of regression coefficients, and *lambda* is the linear predictor (vector of length *n*).

On the regression coefficients we place conditionally Gaussian prior distributions: we assume that, conditionally on the mean vector *b* and covariance matrix *B*,
<p align="center">
(β | b, B) ~ N(b,B).
</p>
With conditional Gaussian prior, we refer to a possibly hierarchical prior, where *b* and/or *B* are random and assigned a hyper-prior. 
Examples include straightforward Gaussian prior distributions with informative(b,B) fixed using prior information, or scale mixtures of Gaussians wherebis set to zeroand the variance has a suitable hierarchical representation such as the Bayesian lasso prior(Park and Casella, 2008), the horseshoe prior, and its extensions (Carvalho et al., 2010; Piironen and Vehtari, 2017).
