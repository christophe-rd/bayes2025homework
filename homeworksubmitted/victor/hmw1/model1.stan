functions {
  real gaussian(real x, real a, real b, real c) {
    return 0.7*exp(-(a*x-b)^2/(2*c^2));
  }
}

data {
  int<lower=1> N;
  vector[N] gdd;
  vector[N] y;
}

parameters {
  real<lower=0> a; // effect size of GDD..?
  real<lower=0> b; // position of the optimum (peak)
  real<lower=0> c; // width of the bell
  real<lower=0> sigma;  
}

transformed parameters {
  real<lower=0> yhat[N];
  
  for(i in 1:N){
    yhat[i] = gaussian(gdd[i],a,b,c);
  }
}

model {
  // Priors
  a ~ normal(1, 0.5);
  b ~ normal(2000, 250);
  c ~ normal(500, 130);
  sigma ~ normal(0, 0.2);
  
  // Model
  y ~ normal(yhat, sigma);
  
}
