#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include <RcppArmadillo.h>
#include <boost/math/special_functions/lambert_w.hpp>

double Mahalanobis(arma::vec x, arma::vec center, arma::mat cov);

const double log2pi = std::log(2.0 * M_PI);

double dmvnorm_arma(arma::vec x, arma::vec mean, arma::mat sigma, bool log);

double logfactorial(int M);

double loglik(arma::vec y, arma::mat X, arma::vec beta);

double logprior(arma::vec beta, arma::vec b, arma::mat B);

double logpost(arma::vec y, arma::mat X, arma::vec beta, arma::vec b, arma::mat B);

double pgmean(double b, double c);

double r_root(double lambda, double d);

#endif


