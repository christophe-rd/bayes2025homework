
data{
  
  int n;
  real<lower=0, upper=1> y[n];
  real<lower=0, upper=1> x[n];
  
}

parameters{
  
  real<lower=0, upper=1> a;
  real<lower=0> d;
  real<lower=0, upper=1> m;
  real<lower=0> phi;
  simplex[2] lambda;
}

model{
  
  vector[n] mu;
  row_vector[2] shape ;
  
  a ~ normal(0.5,0.25);
  d ~ normal(10,5);
  m ~ normal(0.5,0.25);
  phi ~ normal(30,10);
  lambda ~ dirichlet(rep_vector(1,2));
  
  for(i in 1:n){
    
    if(y[i] == 1){
      target += log(lambda[1]);
    } else {
      mu[i] = a/(1+exp(d*(m-x[i])/a));
      shape = [mu[i] * phi, (1 - mu[i] ) * phi];
      target += log(lambda[2]) + beta_lpdf(y[i] | shape[1], shape[2]);
    } 
    
  }
  
}

generated quantities {
  
  real ysim[12];
  real xsim[12];
  
  for (i in 1:10) {
    xsim[i] = uniform_rng(0, 1);
    if(binomial_rng(1, lambda[1]) == 1){
      ysim[i] = 1;
    }else{
      real mu = a/(1+exp(d*(m-xsim[i])/a));
      ysim[i] = beta_rng(mu * phi, (1 - mu) * phi);
    }
  } 
  
  // add some zeros
  for (i in 11:12) {
    xsim[i] = 0;
    if(binomial_rng(1, lambda[1]) == 1){
      ysim[i] = 1;
    }else{
      real mu = a/(1+exp(d*(m-xsim[i])/a));
      ysim[i] = beta_rng(mu * phi, (1 - mu) * phi);
    }
  } 
  
}
