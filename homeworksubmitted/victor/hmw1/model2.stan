functions {
  real gaussian2(real x, real g,real a, real b, real db, real c, real d) {
    return 0.7*exp(-(a*x-(b+g*db)+x*g*d)^2/(2*c^2));
  }
}

data {
  int<lower=1> N;
  vector[N] gdd;
  int<lower=0,upper=1> g[N];
  vector[N] y;
}

parameters {
  real<lower=0> a; // effect size of GDD..?
  real<lower=0> b; // position of the optimum (peak) of population 0
  real db; // position of the optimum (peak) of population 1, relative to pop 0
  real<lower=0> c;  // "width of the bell"
  real d;  // interaction term
  real<lower=0> sigma;
}

transformed parameters {
  real<lower=0> yhat[N];
  // intercept only model
  for(i in 1:N){
    yhat[i] = gaussian2(gdd[i],g[i],a, b,db,c,d);
  }
}

model {
  // Priors
  a ~ normal(1, 0.5);
  b ~ normal(2000, 250);
  db ~ normal(0, 200);
  c ~ normal(500, 130);
  d ~ normal(0, 0.3);
  sigma ~ normal(0, 0.2);
  
  // Model
  y ~ normal(yhat, sigma);
  
}
