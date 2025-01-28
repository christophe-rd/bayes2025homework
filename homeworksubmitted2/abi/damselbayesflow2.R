#Abi's homework 2 for Bayesian modelling January 2025

options(stringsasfactors = FALSE)
library(rstan)
library(rstanarm)
library(dplyr)
library(shinystan)
library(ggplot2)
#I am using the damselfly swimming speed data

fastswim <- read.delim("/bayes2025homework/analyses/input/damselflies/fastswim.dta", sep = "")

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

#My real data does have more species with lower sample sizes, but I want more samples for test data
n.per.sp <- round(runif(nspecies, min = 10, max = 100))

#number of samples total
n <- sum(n.per.sp)

#species level parameters (see prior predictive checks)
mu_int <- 10
sig_int <- 4
intspp <- rnorm(nspecies, mu_int, sig_int)

#while slopes may vary between species, I am not expecting that there is a huge effect of mass nor will it vary a ton between species
mu_slope <- 0.8
sig_slope <- 0.4
slopespp <- rnorm(nspecies, mu_slope, sig_slope)

#overall error
sig_y <- 7
paramsgiven <- c(mu_int, sig_int, mu_slope, sig_slope, sig_y)
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
  sim.data$wetmass[sim.data$species == sp] <- abs(rnorm(n = n.per.sp[sp],m,d))
}
#For simulating my y data, I am getting some nonsensical negative speeds, 
for (o in 1:n){
  sp <- sim.data$species[o]
  sim.data$ypred[o] <- intspp[sp] + slopespp[sp]*sim.data$wetmass[o]
  sim.data$y[o] <- rnorm(n=1, sim.data$ypred[o], sig_y)
}

#compare my real and simulated data
par(mfrow = c(2,1))
plot(fastswim$wetmass, fastswim$speed, main = "real data")
plot(sim.data$wetmass, sim.data$y, main = "simulated data")

# Running the model on simulated data and assessing -----------------------

#simfit_rstanarm <- stan_lmer(y ~ wetmass|species, data = sim.data)

#summary(simfit_rstanarm)
#launch_shinystan(simfit_rstanarm)
#The rstanarm model is giving me a low ESS warning, none of them are actually below 10%, but the intercept is quite close (411 on 1st run)
#my other diagnostics look quite good,all Rhats are 1.0 and no divergent transitions, however, I wanted to actually run my model in Stan to see if I could get better ESS
#The estimates were also still fairly good for my parameters, especially for mean alpha and overall sigma, but it did struggle with the beta means

##stan model edited from example code
#parameters given in stan model are similar to those used to generate data above

##format data
mass <- sim.data$wetmass
species <- as.numeric(as.factor(sim.data$species))
y <- sim.data$y
simfit_stan <- stan("Documents/Github/bayes2025homework/homeworksubmitted2/abi/hierslopeint_editedabi.stan", 
              data=c("n", "nspecies", "y", "mass", "species"), 
              iter=4000, chains=4, seed=377)
launch_shinystan(simfit_stan)
#My ESS was much better for the model in Stan, 
#the rest of my diagnostics are also looking good with well mixed chains, all Rhats below 1.01 (yay)
#Extract estimated parameters from my model and plot
sim_sumar <- summary(simfit_stan)$summary
aparams <- sim_sumar[grep("a|mu_int", rownames(sim_sumar)), c("mean", "2.5%", "25%", "50%", "75%", "97.5%")] %>% as.data.frame()
aparams$names <- rownames(aparams)
bparams <- sim_sumar[grep("b|mu_slope", rownames(sim_sumar)), c("mean", "2.5%","25%", "50%", "75%", "97.5%")] %>% as.data.frame()
bparams$names <- rownames(bparams)

#mean intercept estiamtes for each species and total
ggplot(aparams) +
  geom_point(mapping = aes(x = names, y = mean)) +
  geom_abline(slope = 0, intercept = mu_int, color = "red") +
  labs(title = "intercept estimate means (total and per species), true total mean = 10",
       x = "species groups, mu_int = total mean estimate")
#mean slope estimates for each species and total
ggplot(bparams) +
  geom_point(mapping = aes(x = names, y = mean)) +
  geom_abline(slope = 0, intercept = mu_slope, color = "red") +
  labs(title = "slope estimate means (total and per species), true total mean = 0.9",
       x = "species groups, mu_slope = total mean estimate")

paramsout <- sim_sumar[c("mu_int", "sig_int", "mu_slope", "sig_slope", "sig_y"),"mean"]

#Looking at this comparison table, model was generally able to return the given parameters
paramscompare <- data.frame(given = paramsgiven, estimates = paramsout, row.names = c("mu_int", "sig_int", "mu_slope", "sig_slope", "sig_y"))
print(paramscompare)

# Running the model on empirical data -------------------------------------

#I'm not even going to try to run the rstanarm version on my real data since I don't think ESS problem is going to get better with even less data
#also am running out of time so I did not deal with the fact that I have NAs in my data properly, but I am going to come back to this
fastswim <- fastswim[!is.na(fastswim$speed),]
species <- as.numeric(as.factor(fastswim$species))
mass <- fastswim$wetmass
y <- fastswim$speed
n <- nrow(fastswim)
fit_stan <- stan("Documents/Github/bayes2025homework/homeworksubmitted2/abi/hierslopeint_editedabi.stan", 
                 data=c("n", "nspecies", "y", "mass", "species"), 
                 iter=4000, chains=4, seed=377)
launch_shinystan(fit_stan)
summary(fit_stan)
#My Neffs are all in the thousands, and all Rhat are below 1.01 so I believe everything is working ok! But lets look at how the parameter estimates went
sumar <- summary(simfit_stan)$summary %>% as.data.frame()
paramsout_emp <- sumar[c("mu_int", "sig_int", "mu_slope", "sig_slope", "sig_y"),"mean"]
paramscompare_emp <- data.frame(given = paramsgiven, estimates = paramsout, row.names = c("mu_int", "sig_int", "mu_slope", "sig_slope", "sig_y"))
print(paramscompare_emp)

#Lets plot the predicted fits for each species
plot(1, type="n", ylim=c(-2, 50), xlim=c(0, 60),
     xlab="mass", ylab="swim speed")
for(i in 1:nspecies){
  arow <- i + 5
  brow <- i + 21
  abline(a=sumar$mean[arow], sumar$mean[brow])
}

#Looks pretty close to what I predicted, most species have positive relationship between mass and speed, but there are some species for which it is less important
# Prior predictive checks -------------------------------------------------
#Here I go through my process for prior predictive checks for my parameters, I have given them slightly different names here so as not to mess with my model/data above
#Lets look at our real data
par(mfrow = c(1,1))
plot(fastswim$wetmass, fastswim$speed, colour = fastswim$species)
#slope seems pretty close to 1, and I'll guess that the mean intercept may be somewhere around 15 (and is probably fairly variable)
mu_a <- 15
sig_a <- 4
aspp <- rnorm(nspecies, mu_a, sig_a)

#I'm going to keep the error in slopes fairly low for now
mu_b <- 1
sig_b <- 0.4
bspp <- rnorm(nspecies, mu_b, sig_b)

#overall error
sigma_y <- 7

#Lets make some data, I want to have structure within my data similar to the real data, so I am setting unique means and sd for mass within groups
for(sp in 1:nspecies){
  m <- mass.means$x[sp]
  d <- mass.sd$x[sp]
  sim.data$wetmass[sim.data$species == sp] <- abs(rnorm(n = n.per.sp[sp],m,d))
}
#For simulating my y data, I am getting some nonsensical negative speeds, but I am going to ignore this for now, as I was struggling to fit a gamma distribution 
for (o in 1:n){
  sp <- sim.data$species[o]
  sim.data$ypred[o] <- aspp[sp] + bspp[sp]*sim.data$wetmass[o]
  sim.data$y[o] <- rnorm(n=1, sim.data$ypred[o], sig_y)
}

#compare my real and simulated data
par(mfrow = c(2,1))
plot(fastswim$wetmass, fastswim$speed, main = "real data")
plot(sim.data$wetmass, sim.data$y, main = "simulated data")

#Some of my predicted speeds are quite high compared with the real data, so i'll decrease the slope mean somewhat
mu_b2 <- 0.8
sig_b2 <- 0.4
bspp2 <- rnorm(nspecies, mu_b2, sig_b2)

#redo data simulation
for(sp in 1:nspecies){
  m <- mass.means$x[sp]
  d <- mass.sd$x[sp]
  sim.data$wetmass[sim.data$species == sp] <- abs(rnorm(n = n.per.sp[sp],m,d))
}
for (o in 1:n){
  sp <- sim.data$species[o]
  sim.data$ypred[o] <- aspp[sp] + bspp2[sp]*sim.data$wetmass[o]
  sim.data$y[o] <- rnorm(n=1, sim.data$ypred[o], sig_y)
}

#compare my real and simulated data again
par(mfrow = c(2,1))
plot(fastswim$speed~fastswim$wetmass, main = "real data")
plot(sim.data$y~sim.data$wetmass, main = "simulated data")

#the relationship is looking better, but I think my average speed is looking a little high, (especially for the smaller individuals) so I'm going to lower the mean intercept
mu_a2 <- 10
sig_a2 <- 4
aspp2 <- rnorm(nspecies, mu_a2, sig_a2)

#One more time!!
for(sp in 1:nspecies){
  m <- mass.means$x[sp]
  d <- mass.sd$x[sp]
  sim.data$wetmass[sim.data$species == sp] <- abs(rnorm(n = n.per.sp[sp],m,d))
}
for (o in 1:n){
  sp <- sim.data$species[o]
  sim.data$ypred[o] <- aspp2[sp] + bspp2[sp]*sim.data$wetmass[o]
  sim.data$y[o] <- rnorm(n=1, sim.data$ypred[o], sig_y)
}

#compare my real and simulated data again
par(mfrow = c(2,1))
plot(fastswim$speed~fastswim$wetmass, main = "real data")
plot(sim.data$y~sim.data$wetmass, main = "simulated data")

#I think I'm happy with these
#Histogram of possible intercepts (adding variation to estimates)
hist(rnorm(10000, mu_a2, sig_a2))

ppcreps <- 200
aprior <- rnorm(ppcreps, mu_a2, sig_a2)

plot(1, type="n", ylim=c(-2, 50), xlim=c(0, 60),
     xlab="mass", ylab="swim speed")
for (i in c(1:ppcreps)){
  abline(a=aprior[i], b=mu_b2)
}

#this looks ok, I am getting a small number of negative intercepts though which is not great, this parameter should probably be fit with a different distribution (gamma or lognormal perhaps)
hist(rnorm(1000, mu_b2, sig_b2))
bprior <- rnorm(ppcreps, mu_b2, sig_b2)

plot(1, type="n", ylim=c(-2, 50), xlim=c(0, 60),
     xlab="mass", ylab="swim speed")
for (i in c(1:ppcreps)){
  abline(a=mu_a2, b=bprior[i])
}
#This looks fine I think! I'm generally expecting my slopes to be slightly positive, but I also do not think its unreasonable to believe that for some species they could get slower as they get bigger, 
#However, negative slopes would predict negative speeds at higher masses, but this is more a limitation of using a linear model (since it also is not realistic for species speed to increase infinitely with increasing mass)
#And very few of my slopes are negative so I'm not super concerned for the purposes of this assignment
