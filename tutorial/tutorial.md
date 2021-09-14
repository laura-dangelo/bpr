## `bpr` tutorial

This tutorial is a short guide on the use of the `bpr` package to perform posterior sampling and inference on the parameters of Bayesian Poisson regression models.

### Model specification
First, we briefly recall the model specification and prior distributions available in the package.
The general setting is that of standard Poisson regression models, where the interest is to regress a vector of counts *y*, of size *n*, on an *n x p* matrix *X* of covariates:

y ~ Poisson(lambda); lambda = exp{ X * beta }

where *beta* is the *p x 1* vector of regression coefficients.
