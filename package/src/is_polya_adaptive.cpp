#include "commonfunctions.h"

// [[Rcpp::export(.ISPolya)]]
Rcpp::List ISPolya(const int n_rep, 
                    arma::vec y, arma::mat X,
                    arma::vec b, arma::mat B,
                    double r_start, 
                    arma::vec beta_start,
                    double max_dist,
                    int max_r,
                    double trunc_lambda) 
{
  int p = X.n_cols ;
  int n = X.n_rows ;
   
  // IS quantities 
  arma::vec w(n_rep) ;
  arma::vec logw(n_rep) ;
  
  arma::vec logr(n) ;
  
  arma::vec dif = (y-r_start)/2 ;
  arma::mat invB = arma::inv_sympd(B) ;
  
  arma::mat beta_out = arma::zeros(p, n_rep) ;
  beta_out.col(0) = beta_start ;
  arma::mat r_out = arma::zeros(n, n_rep) ;
  r_out.col(0).fill(r_start) ;
  
  arma::vec beta = beta_start ;
  arma::vec beta_center = beta_start ;
  
  arma::vec r_new(n) ;
  int r_tmp ;
  arma::vec lambda_truncated(n) ;
  
  // mean and variance (conditioned on old beta)
  arma::vec lin_pred(n) ;
  arma::vec omega_mean(n) ;
  arma::mat V(n,n) ;
  arma::vec m(n) ;
  arma::vec kappa(n) ;
  

  for(int i = 0; i < n_rep ; i++)
  {
    
    lin_pred = X * beta_center ; // linear predictor
    for(int j = 0; j < n; j++)
    {
      lambda_truncated(j) = std::min( exp(lin_pred(j)), trunc_lambda ) ;
      r_tmp = ceil( r_root(lambda_truncated(j), max_dist) ) ;
      r_new(j) = std::min(r_tmp, max_r) ;
    }
    logr = log( r_new ) ;
    
    for(int j = 0; j < n; j++)
    {
      omega_mean(j) = pgmean(r_new(j) + y(j), lin_pred(j) - logr(j)) ;
      dif(j) = ( y(j) - r_new(j) )/2 ;
    }
    
    V = arma::inv_sympd( X.t() * arma::diagmat(omega_mean) * X + invB ) ;
    for(int j = 0; j < n; j++)
    {
      kappa(j) = omega_mean(j) * logr(j) + dif(j) ;
    }
    m = V * ( X.t() * kappa + invB * b ) ;
    
    // new beta
    beta = arma::mvnrnd(m, V) ; 
    
    if( logpost(y, X, beta, b, B) - logpost(y, X, beta_center, b, B) > 0 ) { beta_center = beta ; }
    
    // weights
    w(i) = logpost(y, X, beta, b, B) - dmvnorm_arma(beta, m, V, true) ;
    logw(i) = logpost(y, X, beta, b, B) - dmvnorm_arma(beta, m, V, true) ;
    
    beta_out.col(i) = beta ; 
    r_out.col(i) = r_new ;
    
  }
    
  for(int t = 0; t < n_rep; t++)
  {
    w(t) = exp( w(t) ) ;
  }
  
  return Rcpp::List::create(Rcpp::Named("beta") = beta_out.t(),
                            Rcpp::Named("r") = r_out.t(),
                            Rcpp::Named("logw") = logw,
                            Rcpp::Named("w") = w
                           );
}
