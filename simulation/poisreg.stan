data {
  int<lower=1> N;
  int<lower=1> p;   // number of predictors
  matrix[N, p] X;   // predictor matrix
  int<lower=0> y[N];      // outcome vector
}
parameters {
  vector[p] beta;       // coefficients for predictors
}
model {
  for(j in 1:p) beta[j] ~ normal(0, 2);
  y ~ poisson_log(X * beta);
}