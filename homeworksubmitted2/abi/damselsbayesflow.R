#Abi's homework 2 for Bayesian modelling January 2025

options(stringsasfactors = FALSE)
library(rstan)
library(rstanarm)
library(dplyr)
library()
#I am using the damselfly swimming speed data
fastswim <- read.delim("bayes2025homework/analyses/input/damselflies/fastswim.dta", sep = "")
#for some reason, some variables read in as characters so I need to change the classes
for( c in 4:13) {
  fastswim[,c] <- as.numeric(fastswim[,c])
}

# Model Building ----------------------------------------------------------

#I am predicting that there is going to be a generally linear relationship between swim speed and body mass
#I also am predicting that both the minimum swim speed (intercept) and the influence of mass on speed (slope) will vary between species
#Mass will be positively linearly related to wet mass
#It does not make sense for either mass or swim speed to be negative, and both of these values are close to 0 in the data, so I will not be able to fit a normal distribution, I went with gamma instead.
miu.alpha <- 7
sig.alpha <- 4
alphaspp <- rnorm(miu.alpha, sig.alpha)

shape.beta <- 0.7
scale.beta <- 0.3
betaspp <- rnorm(miu.beta, sig.beta)


speedfit <- stan_glm(fastswim$speed ~ alphaspp + betaspp*fastswim$mass, data = fastswim)
# Simulating data ---------------------------------------------------------

#I want my simulated data to have the same number of grouping factors (species) as my real data
nspecies <- length(unique(fastswim$species))

#species level parameters, setting the distributions for the intercepts and slopes of the model
miu.alpha <- 15
sig.alpha <- 7
alphaspp <- rnorm(nspecies, miu.alpha, sig.alpha)

miu.beta <- 0.7
sig.beta <- 0.3
betaspp <- rnorm(nspecies, miu.beta, sig.beta)

sig.y <- 10
#creating a similar-ish distribution of sample sizes per species (more species with less samples)
real.n.per.sp <- dplyr::count(fastswim, fastswim$species)
n.per.sp <- round(rpois(nspecies, mean(real.n.per.sp$n)))
#number of samples total
n <- sum(n.per.sp)
#create data frame
sim.data <- data.frame(species = rep(1:nspecies, n.per.sp), 
                       mass = NA,
                       speed = NA)
#simulation group means and var, based on real max and mins of groups
mass.means <- aggregate(fastswim$wetmass, list(fastswim$species), FUN = "mean")
mass.sd <- aggregate(fastswim$wetmass, list(fastswim$species), FUN = "sd")
#For generating masses, I used the absolute value of my draws from rnorm, 
#I realize this is increasing the means and decreasing the variance of my distributions somewhat, but it was a convenient way to prevent generating negative masses.
  for(sp in 1:nspecies){
  m <- mass.means$x[sp]
  d <- mass.sd$x[sp]
  sim.data$mass[sim.data$species == sp] <- abs(rnorm(n = n.per.sp[sp],m,d))
}
  for (o in 1:n){
    sp <- sim.data$species[o]
   yhat <- alphaspp[sp] + betaspp[sp]*sim.data$mass[o]
    sim.data$speed[o] <- rnorm(1, yhat, sig.y)
  }
plot(fastswim$wetmass, fastswim$speed)
plot(sim.data$mass, sim.data$speed)
# Running the model on simulated data and assessing -----------------------

speedfit <- stan_glm(sim.data$speed ~ alphaspp + betaspp*sim.data$mass,)
