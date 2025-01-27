#Abi's homework 2 for Bayesian modelling January 2025

options(stringsasfactors = FALSE)
library(rstan)
library(rstanarm)
library(dplyr)
library(shinystan)
#I am using the damselfly swimming speed data

fastswim <- read.delim("~/Documents/GitHub/bayes2025homework/analyses/input/damselflies/fastswim.dta", sep = "")

#for some reason, some variables read in as characters so I need to change the classes
for( c in 4:13) {
  fastswim[,c] <- as.numeric(fastswim[,c])
}

# Model Building ----------------------------------------------------------

#I am predicting that there is going to be a generally linear relationship between swim speed and body mass
#I also am predicting that both the minimum swim speed (intercept) and the influence of mass on speed (slope) will vary between species
#Mass will be positively linearly related to wet mass

#I want my simulated data to have the same number of grouping factors (species) as my real data
nspecies <- length(unique(fastswim$species))

#My real data does have more species with lower sample sizes, but I will ignore this for now
n.per.sp <- round(runif(nspecies, min = 10, max = 50))
hist(n.per.sp)

#number of samples total
n <- sum(n.per.sp)

mu.alpha <- 10
sig.alpha <- 3.5

intspp <- rnorm(n, mu.alpha, sig.alpha)

mu.beta <- 1.5
sig.beta <- 0.9
slopespp <- rnorm(n, mu.beta, sig.beta)

#overall error
sig.y <- 5

# Simulating data ---------------------------------------------------------

#create data frame
sim.data <- data.frame(species = rep(1:nspecies, n.per.sp))
sim.data$species <- as.factor(sim.data$species)
#simulation group means and var, based on real max and mins of groups, I wanted each group to have its own internal distribution as this would be most similar to my actual data structure
mass.means <- aggregate(fastswim$wetmass, list(fastswim$species), FUN = "mean")
mass.sd <- aggregate(fastswim$wetmass, list(fastswim$species), FUN = "sd")

#For generating masses, I used the absolute value of my draws from rnorm
#I realize this is increasing the means and decreasing the variance of my distributions somewhat, but it was a convenient way to prevent generating negative masses.
for(sp in 1:nspecies){
  m <- mass.means$x[sp]
  d <- mass.sd$x[sp]
  sim.data$mass[sim.data$species == sp] <- abs(rnorm(n = n.per.sp[sp],m,d))
}
for (o in 1:n){
  sp <- sim.data$species[o]
  sim.data$ypred[o] <- intspp[sp] + slopespp[sp]*sim.data$mass[o]
  sim.data$y[o] <- rnorm(1, sim.data$ypred[o], sig.y)
}

#compare my real and simulated data
par(mfrow = c(2,1))
plot(fastswim$wetmass, fastswim$speed, main = "real data")
plot(sim.data$mass, sim.data$y, main = "simulated data")

# Running the model on simulated data and assessing -----------------------

testfit <- stan_lmer(y ~ mass|species,
                       data = sim.data)

summary(speedfit)

launch_shinystan(speedfit)


# Running the model on empirical data -------------------------------------

realfit <- stan_lmer(speed ~ wetmass|species,
                     data = fastswim,
                     iter = 5000)
