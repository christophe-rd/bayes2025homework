//
// This Stan program defines a simple model, with a
// vector of values 'yobs' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=0> N;       // Number of observations
  vector[N] x;          // Predictor
  vector[N] y;          // Response
}
parameters {
  real a;           // Intercept
  real beta1;           // Slope
  real<lower=0> sigma;  // Error SD
}
model {
  y ~ normal(a + beta1 * x, sigma); // Likelihood
}



