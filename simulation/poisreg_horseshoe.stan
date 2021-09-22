data {
  int<lower=1> N;
  int<lower=1> p;   // number of predictors
  matrix[N, p] X;   // predictor matrix
  int<lower=0> y[N];      // outcome vector
  real tau;
}
parameters {
  vector[p] beta;       // coefficients for predictors
  vector<lower=0>[p] lambda;
  // real<lower=0> tau;
}
model {
  for(j in 1:p) lambda[j] ~ cauchy(0, 1);
  // tau ~ cauchy(0, 1) ;
  
  for(j in 1:p) beta[j] ~ normal(0, tau * fabs(lambda[j]) );
  y ~ poisson_log(X * beta);
}
