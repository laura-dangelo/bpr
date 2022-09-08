# `bpr` package
## Posterior sampling and inference for Bayesian Poisson Regression

Efficient C++ based R package to sample from the posterior distribution of Poisson regression models. 
The model specification makes use of Gaussian (or conditionally Gaussian) prior distributions on the regression coefficients. 

A tutorial on how to use the package can be found [here](/tutorial).

### Installation
The package is available on CRAN. To install it, simply run
```r
install.packages("bpr")
```

#### Alternative installation from Github
```r
library(devtools)
install_github("laura-dangelo/bpr/package")
```

The package requires [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html#:~:text=The%20'Rcpp'%20package%20provides%20R,integration%20of%20third%2Dparty%20libraries.), [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) and the C++ library [boost](https://www.boost.org/).



## How to cite
Thank you for your interest in my work! If you use this package in any of your projects, please cite the package and the related paper as:

- D’Angelo, L. (2021), 'bpr: Bayesian Poisson regression', URL: https://CRAN.R-project.org/package=bpr

- D'Angelo, L. and Canale, A. (2022), 'Efficient posterior sampling for Bayesian Poisson regression', arXiv preprint at _arXiv:2109.09520_
