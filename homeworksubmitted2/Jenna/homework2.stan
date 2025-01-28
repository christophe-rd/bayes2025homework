//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; 	// number of total observations
  int<lower=0> Nspp; 	// number of species (grouping factor)
  int sppID[N]; 	// species identity, coded as int
  vector[N] centered_latitude; 	// latitude of data point (predictor for slope)
  real CsupL[N]; 		// tooth length (response)
}

// The parameters accepted by the model.
parameters {
  real mu_a;		// mean intercept across species
real<lower=0> sigma_a;	// variation of intercept across species	
real mu_b;		// mean slope across species
real<lower=0> sigma_b;	// variation of slope across species
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
real Aspp[Nspp]; 		//the intercept for each species
real Bspp[Nspp]; 		//the slope for each species 
}


transformed parameters{
real ypred[N];
for (i in 1:N){
    ypred[i]=Aspp[sppID[i]]+Bspp[sppID[i]]*centered_latitude[i];
}
}


model{	
Bspp ~ normal(mu_b, sigma_b); // this creates the partial pooling on slopes
Aspp ~ normal(mu_a, sigma_a); // this creates the partial pooling on intercepts
CsupL ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
// Priors ...
mu_a ~ normal(5,5);
sigma_a ~ normal(0,5);
mu_b ~ normal(0,1);
sigma_b ~ normal(0,1);
sigma_y ~ normal(0,1);
}	

generated quantities {
  vector[N] y_rep;        // Posterior predictive samples
  for (n in 1:N) {
    y_rep[n] = normal_rng(ypred[n], sigma_y);  // Simulate new data points
  }
}

