#include "commonfunctions.h"

// [[Rcpp::export(.MHPolya_adaptr)]]
Rcpp::List MHPolya_adaptr(const int n_rep, 
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
  
  // MH quantities
  int accepted = 0 ;
  double acc ;
  double alpha ;
  double u ;
  
  arma::vec logr(n) ;
  
  arma::vec dif = (y-r_start)/2 ;
  arma::mat invB = arma::inv_sympd(B) ;
  
  arma::mat beta_out = arma::zeros(p, n_rep) ;
  beta_out.col(0) = beta_start ;
  arma::mat r_out = arma::zeros(n, n_rep) ;
  r_out.col(0).fill(r_start) ;
  
  arma::vec beta_old = beta_out.col(0) ;
  arma::vec beta_new(p) ;
  
  arma::vec r_new(n) ;
  int r_tmp ;
  arma::vec lambda_truncated(n) ;
  
  // media e varianza condizionate al vecchio beta
  arma::vec lin_pred(n) ;
  arma::vec omega_mean(n) ;
  arma::mat V(n,n) ;
  arma::vec m(n) ;
  arma::vec kappa(n) ;
  
  // media e varianza condizionate al nuovo beta
  arma::vec lin_pred_new(n) ;
  arma::vec omega_mean_new(n) ;
  arma::mat V_new(n,n) ;
  arma::vec m_new(n) ;
  arma::vec kappa_new(n) ;
  
  
  for(int i = 0; i < n_rep - 1; i++)
  {
    
    lin_pred = X * beta_old ; // linear predictor
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

    
    // estraggo nuovo valore per beta
    beta_new = arma::mvnrnd(m, V) ; 
    lin_pred_new = X * beta_new ;

    for(int j = 0; j < n; j++)
    {
      omega_mean_new(j) = pgmean(r_new(j) + y(j), lin_pred_new(j) - logr(j)) ;
    }
    V_new = arma::inv_sympd(X.t() * arma::diagmat(omega_mean_new) * X + invB) ;
    
    for(int j = 0; j < n; j++)
    {
      kappa_new(j) = omega_mean_new(j) * logr(j) + dif(j) ;
    }
    m_new = V_new * (X.t() * kappa_new + invB * b) ;

    
    alpha = exp( logpost(y, X, beta_new, b, B) - 
      logpost(y, X, beta_old, b, B) +
      dmvnorm_arma(beta_old, m_new, V_new, true) -
      dmvnorm_arma(beta_new, m, V, true) ) ;
    
    u = ::Rf_runif(0, 1) ;
    if(u < alpha)
    {
      beta_old = beta_new ;
      accepted++ ;
    }
    beta_out.col(i+1) = beta_old ; 
    r_out.col(i+1) = r_new ;
    acc = (1. * accepted) / n_rep ;

  }
  return Rcpp::List::create(Rcpp::Named("beta") = beta_out.t(),
                            Rcpp::Named("r") = r_out.t(),
                            Rcpp::Named("acceptance_rate") = acc);
}
