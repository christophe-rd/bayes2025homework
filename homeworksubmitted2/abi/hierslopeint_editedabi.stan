// Two-level (1 hierarchical grouping) `random' slope and intercept model
// Partial pooling on intercepts and slopes 

data{
int<lower=0> n; 	// number of total observations
int<lower=0> nspecies; 	// number of species (grouping factor)
int species[n]; 	// species identity, coded as int
vector[n] mass; 	// wetmass of individual (predictor for slope)
real y[n]; 		// max swim speed of individual (avg of 3) (response)
}

parameters{
real mu_int;		// mean intercept across species
real<lower=0> sig_int;	// variation of intercept across species	
real mu_slope;		// mean slope across species
real<lower=0> sig_slope;	// variation of slope across species
real<lower=0> sig_y; 	// measurement error, noise etc. 	
real a[nspecies]; 		//the intercept for each species
real b[nspecies]; 		//the slope for each species 

}

transformed parameters{
real ypred[n];
for (i in 1:n){
    ypred[i]=a[species[i]]+b[species[i]]*mass[i];
}
}

model{	
b ~ normal(mu_slope, sig_slope); // this creates the partial pooling on slopes
a ~ normal(mu_int, sig_int); // this creates the partial pooling on intercepts
y ~ normal(ypred, sig_y); // this creates an error model where error is normally distributed
// Priors ...
mu_int ~ normal(10,2);
sig_int ~ normal(0,4);
mu_slope ~ normal(0.8,0.2);
sig_slope ~ normal(0,0.4);
sig_y ~ normal(0,7);
}	
