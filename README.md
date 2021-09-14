# BPR package
## Posterior sampling and inference for Bayesian Poisson Regression

Efficient C++ based R package to sample from the posterior distribution of Poisson regression models. 
The model specification makes use of Gaussian (or conditionally Gaussian) prior distributions on the regression coefficients. 


### Installation
#### From Github
```{R}
library(devtools)
install_github("laura-dangelo/bpr/package")
```

#### From source
Download the tarball archive "bpr_1.0.tar.gz". In R, run
```{R}
install.packages("path/to/bpr_1.0.tar.gz", repos = NULL, source = TRUE)
```

The package requires [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html#:~:text=The%20'Rcpp'%20package%20provides%20R,integration%20of%20third%2Dparty%20libraries.), [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) and the C++ library [boost](https://www.boost.org/).
