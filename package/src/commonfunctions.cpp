#include "commonfunctions.h"

using boost::math::lambert_w0;

// Multivariate normal density
double Mahalanobis(arma::vec x, arma::vec center, arma::mat cov) {
  arma::vec x_cen;
  x_cen.copy_size(x);
  x_cen = x - center;
  return arma::accu((x_cen.t() * cov.i()) % x_cen.t());    
}

double dmvnorm_arma(arma::vec x, arma::vec mean, arma::mat sigma, bool log = false) { 
  double distval = Mahalanobis(x, mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  double logretval = -( (x.n_rows * log2pi + logdet + distval)/2  ) ;
  if (log) { 
    return(logretval);
  } else { 
    return(exp(logretval));
  }
} 

// factorial function
double logfactorial(int M){
  double fact = 0;
  if((M == 1) | (M == 0)) return 0;
  for(int l = 1; l <= M; l++ ){
    fact += std::log(l);
  }
  return fact;
}

// Poisson log likelihood
double loglik(arma::vec y, arma::mat X, arma::vec beta){
  int n = y.n_elem;
  int k = 0;
  arma::vec loglambda = X * beta;
  arma::vec llik(n);
  
  for(k = 0; k < n; k++) { 
    llik(k) = -std::exp(loglambda[k]) + y[k] * loglambda[k] - logfactorial(y[k]) ; 
  }
  return arma::accu(llik);
}

// log prior
double logprior(arma::vec beta, arma::vec b, arma::mat B){
  return dmvnorm_arma(beta, b, B, true);
}

// log posterior
double logpost(arma::vec y, arma::mat X, arma::vec beta, arma::vec b, arma::mat B){
  return (loglik(y, X, beta) + logprior(beta, b, B));
}

// polya-gamma mean
double pgmean(double b, double c)
{
  return b/(2*c)*(std::exp(c)-1)/(std::exp(c)+1);
}

double r_root(double lambda, double d)
{
  double c = exp( lambda ) * (d + 1) ;
  double logc = log(c) ;
  double r = lambda * logc /  ( logc + lambda * lambert_w0(- pow(c, -1/lambda) * logc /lambda ) ) ;
  return(r) ;
}

