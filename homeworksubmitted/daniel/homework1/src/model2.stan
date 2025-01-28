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
  int<lower=0> N;         // Number of observations
  vector[N] x1;           // Continuous predictor
  int<lower=0,upper=1> x2[N]; // Binary covariate
  vector[N] y;            // Response variable
}
parameters {
  real a;                 // Intercept
  real b1;                // Effect of x1
  real b2;                // Effect of x2
  real b3;                // Interaction effect (x1 * x2)
  real<lower=0> sigma;    // Error SD
}
transformed parameters {
  vector[N] x2_real = to_vector(x2); // Convert integer array to real vector
}
model {
  y ~ normal(a + b1 * x1 + b2 * x2_real + b3 * (x1 .* x2_real), sigma); // Likelihood
}

