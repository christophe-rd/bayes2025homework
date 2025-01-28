### FRST 507 - HOMEWORK 2
## Hannah Bates - Jan 2025

# load pkgs
library(rstanarm)
library(dplyr)
library(tidyr)

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)
par(mfrow=c(1, 1)) # default 1 plot


## IF RUNNING FROM REPO: Set working directory where carnivore data is stored:

# setwd("/temporalecologylab/bayes2025homework/analyses/input/carnivores")
# And update the 2 paths below!


## STEP 1: COME UP WITH A MODEL & prep ----

# I have chosen to work with the carnivore mass & teeth length data. I would 
# like to model how premolar length in carnivores can be predicted by 
# body mass. I would expect a positive linear relationship, since most body parts
# increase proportionally to body size, which increases proportionally to body mass.
# However, this relationship may differ among phylogeny Families: the premolar tooth,
# important for tearing and crushing food, may be more or less important depending
# on the diet of the species, which I am assuming is more similar across species
# in the same Family than between species from different Families.
# Thus, I will group data by Family (1 hierarchical grouping).
# In lmer terms: PM4length ~ mass + (mass | Family). As shown below, both PM4 length
# and body mass become log-transformed. I expect my slopes and intercepts to be
# 'randomly' distributed, rather than fixed.

# PREP
# Load carnivore data
data_mass <- read.csv("INPUT/carnivorebodymass.csv") # UPDATE
data_teeth <- read.csv("INPUT/carnivoreteeth.csv") # UPDATE

# Inspect files
length(unique(data_teeth$Species)) # checking to see how many species
length(unique(data_mass$species))

#Remove NA rows for PM4
data_teeth %>% drop_na(PM4)

# Average the teeth data by species, as we only have averaged mass data
data_teethsum <- data_teeth %>% 
  drop_na(PM4) %>% # remove entries without PM4 length data
  group_by(Family, Species, sex) %>%
  summarize(PM4_length = mean(PM4))

# Clean and combine data based on species and sex
data <- left_join(data_teethsum, data_mass, by = join_by(Species == species, sex)) %>%
  drop_na(mass.midpoint) # remove entries that lack mass data
length(unique(data$Species)) # left with 152 species that have both teeth and mass data
length(unique(data$Family.x)) # left with 11 unique Families... different for mass data for some reason

# Exploratory plotting
cols <- sample(colors(), length(unique(data$Family.x)), replace = FALSE) # colours based on groupings
plot(data$log.midpoint.mass, data$PM4_length, col = cols) # looks like it might need to be transformed
plot(data$log.midpoint.mass, log(data$PM4_length), col = cols) # better!
plot(data_teeth$Long, data_teeth$Lat, col = cols) # A map! Yay for maps!

# Log-transforming PM4 length... not sure if this is the right step to do this...
data$log.PM4 <- log(data$PM4_length)
plot(data$log.midpoint.mass, data$log.PM4, col = cols) # nice!
hist(data$log.PM4)

## STEP 2: SIMULATE TEST DATA ----

# Here, I simulate data similar to my empirical data with set parameters. This way I
# will be able to see if my model is capable of returning parameters similar to 
# the 'actual' values which I will set. This will help me evaluate my model
# performance before trusting any estimates from my empirical data.

# Create family-level parameters
# NOTE: these have been adjusted since my initial simulation based on my PPCs in step 3.
Nfam <- 11 # this is how many I will be working with
mu_mass <- 0.5 # intercept mean
sigma_mass <- 2 # intercept variance
mu_slope <- 0.5 # slope mean
sigma_slope <- 0.3 # I expect the slope to vary considerably between families

family_mass <- rnorm(Nfam, mu_mass, sigma_mass) # intercept hyperprior
family_trend <- rnorm(Nfam, mu_slope, sigma_slope) # slope hyperprior
hist(family_mass)
hist(family_trend)

# Define overall error
sigma_y <- 0.5

# Keeping parameters together
paramsgiven <- c(mu_mass, mu_slope, sigma_slope, sigma_mass, sigma_y)

# Simulate some sweet sweet data
n_data_per_family <- round(runif(Nfam,1,20)) # amount of measurements in each family
family <- rep(1:Nfam, n_data_per_family) # generate random families
N <- length(family) # sample size across all families
mass <- rep(NA, N) # empty vector for loop

# build x data
for (f in 1:Nfam) {
  mass[family==f] <- runif(n_data_per_family[f],2,6) # sample from within the limits of my actual data
}

# build y data
ypred <- length(N)
for (n in 1:N){
  fam <- family[n]
  ypred[n] <- family_mass[fam] + family_trend[fam]*mass[n]
}

yobs <- rnorm(N, ypred, sigma_y) # incorporate overall variance

# Visualize
plot(mass,ypred,xlim = c(0,6), ylim = c(0,4)) # can differentiate between families
plot(mass,yobs,xlim = c(0,6), ylim = c(0,4)) # adjusted sigma_y so that cannot visually differentiate between families
hist(ypred)
hist(yobs) # somewhat normally distributed?

# Test model with simulated data
m1 <- stan_lmer(yobs ~ mass + (mass|family))
estimatedmeans <- summary(m1)[c("(Intercept)","mass","Sigma[family:mass,mass]","Sigma[family:(Intercept),(Intercept)]","sigma"),c("mean","10%","90%")] # store estimates

# Compare with given parameters
paramsgiven # Set
estimatedmeans # Estimated

summary(m1) # Model converged with Rhat = 1 and high enough Neff values

# Plot estimated vs given parameters (not comparable between them but nice to see them all at once)
plot(estimatedmeans[,1], xaxt = "n", 
     main = "Given vs estimated parameters from simulated data",
     xlab = "parameters",
     ylab = "mean parameter estimates", 
     ylim = c(-1,5))
axis(1, at = 1:5, labels = c("mu_int", "mu_slope", "sigma_slope", "sigma_int", "sigma_y"))
points(1:5, paramsgiven, pch = 10)
legend("topright",legend = "Given parameters", pch = 10)
arrows(1:5, estimatedmeans[,2], 1:5, estimatedmeans[,3], angle = 90, code = 3, length = 0.1)

# Results: 
# Estimates: The model managed to return my given estimates fairly well, with all
# given parameters being within the 80% CIs of the estimates except for the mean 
# slope and slope variance. These estimates were precise but their CIs did not
# encompass the given value, so I cannot expect my model to predict these values
# to a certain accuracy. On some runs they do fall within the CIs.
# It was the intercept and variance between groups of the intercept
# which the model struggled to predict precisely, with large error bars.
# Convergence: My model converged, with no errors after running the model. I obtained
# an Rhat value of 1.0 for all estimated parameters, and a Neff >> 10% of the 4000 
# iterations for all parameters, indicating the 4 chains converged.



## STEP 3: PRIOR PREDICTIVE CHECKS ----

# Here, I will make sure the priors I am setting are realistic for what I am trying
# to model.

# These are the initial parameter estimates I gave:
mu_mass <- 0.5 # intercept mean
sigma_mass <- 1 # intercept variance
mu_slope <- 0.5 # slope mean
sigma_slope <- 0.5 # slope variance

# Explore the slope (beta) prior
hist(rnorm(10000, mu_slope, sigma_slope))
# I would expect very few slopes above 2 and below 0, so I will decrease sigma
sigma_slope <- 0.3
hist(rnorm(10000, mu_slope, sigma_slope)) # this looks more similar to what I would expect

ppcreps <- 500 # no of reps
bprior <- rnorm(ppcreps, mu_slope, sigma_slope)
plot(1, type="n", ylim=c(0, 10), xlim=c(0, 10),
     xlab="log(mass)", ylab="log(length)")
for (i in c(1:ppcreps)){ # plot variety of possible slopes
  abline(a=mu_mass, b=bprior[i])
}
# this plots a variety of possible slopes all with the same intercept.
# Overall, these span the expected realm of possibility that I would personally
# expect: that there is usually a positive relationship between mass and premolar length.
# I am not sure if a negative trend would exist, but here I am leaving it as a possibility.

# Explore the intercept (alpha) prior
hist(rnorm(10000, mu_mass, sigma_mass))
# I want to increase the variance as I'm not confident on the boundaries and want to
# allow for more possibilities.
sigma_mass <- 2

aprior <- rnorm(ppcreps, mu_mass, sigma_mass)
plot(1, type="n", ylim=c(0, 7), xlim=c(0, 7),
     main = "PPC: Simulating intercept and slope parameters",
     xlab="log(mass)", ylab="log(length)")
for (i in c(1:ppcreps)){ # plot variety of possible slope AND intercept combos
  abline(a=aprior[i], b=bprior[i])
}
# Overall, these span the expected realm of possibility that I would personally
# expect: that there is usually a positive relationship between mass and premolar length.
# There are some near-flat lines where this relationship may be considerably weak,
# and some high intercepts for Families with generally larger masses.

# I now adjust my parameters in step 2 and make sure the model still works with
# the more explicitly defined priors.

m2 <- stan_lmer(yobs ~ mass + (mass|family),
                prior = normal(mu_slope, sigma_slope), # this time I am defining priors
                prior_intercept = normal(mu_mass, sigma_mass))
summary(m2)
prior_summary(m2)

# Compare with given parameters
estimatedmeans2 <- summary(m2)[c("(Intercept)","mass","Sigma[family:mass,mass]","Sigma[family:(Intercept),(Intercept)]","sigma"),c("mean","10%","90%")] # store estimates
paramsgiven # Set
estimatedmeans2 # Estimated
# Similar results as before


## STEP 4: RUN MODEL ON EMPIRICAL DATA ----

# I will now run the model on my empirical data, using the priors I have defined above.

mreal <- stan_lmer(log.PM4 ~ log.midpoint.mass + (log.midpoint.mass|Family.x),
                   data = data,
                   prior = normal(mu_slope, sigma_slope),
                   prior_intercept = normal(mu_mass, sigma_mass))
summary(mreal)
coef(mreal)

# It looks like the model converged, with Rhat = 1.0 and sufficient Neff (>> 400)
# for all parameters. Note: I've run this a few different times, and occasionally
# I have gotten 1-2 divergent transitions.
# The estimated group-level slopes and intercepts certainly vary,
# but are more similar than what I was expecting. I might go back and reduce my variance
# (sigmas) for slope and intercept based on these findings. To see these differences,
# I have made the plots below.


# Basic plotting with overall trend
plot(data$log.midpoint.mass, data$log.PM4, col = cols,
     main = "Partially-pooled trend of mass vs premolar length",
     xlab = "log(mass)",
     ylab = "log(Pm4 length)")
abline(a = summary(mreal)["(Intercept)","mean"], b = summary(mreal)["log.midpoint.mass","mean"], lwd = 2)

# Plotting lines for each family
plot(data$log.midpoint.mass, data$log.PM4, col = cols, ylim = c(1,4), xlim = c(2,5.5),
     main = "Model-produced trendlines over empirical data",
     xlab = "log(mass)",
     ylab = "log(Pm4 length)")
abline(a = summary(mreal)["(Intercept)","mean"], b = summary(mreal)["log.midpoint.mass","mean"], lwd = 2)
ints <- coef(mreal)$Family.x[1]
slopes <- coef(mreal)$Family.x[2]
for (f in 1:Nfam) {
  abline(a = ints[f,1], b = slopes[f,1], col = cols[f])
}
legend("topleft", legend = unique(data$Family.x), col = cols, pch = 19) # Add a legend


## STEP 5: RETRODICTIVE CHECKS ? ----

# I won't do this step as formally, but here I will see what I can learn from my
# model and output and see if any adjustments to my model should be made.

launch_shinystan(mreal)

# Chains show good mixing, seem to properly explore the parameter space.

hist(posterior_linpred(mreal)) # histogram of posterior distribution - nothing stands out

# Plot estimated slope parameters for each Family
cis <- posterior_interval(mreal, prob = 0.5)[seq(4, 24, by = 2),] + summary(mreal)["log.midpoint.mass","mean"] # 50% CIs
means <- summary(mreal)[seq(4, 24, by = 2),"mean"] + summary(mreal)["log.midpoint.mass","mean"]
plot(means, xaxt = "n", 
     main = "Estimated slope parameters from real data",
     xlab = "Family",
     ylab = "mean slope estimates with 50% CIs", 
     pch = 19,
     ylim = c(0.4,0.8),
     col = cols)
axis(1, at = 1:Nfam, labels = unique(data$Family.x))
arrows(1:Nfam, cis[,1], 1:Nfam, cis[,2], angle = 90, code = 3, length = 0.1)
abline(a=summary(mreal)["log.midpoint.mass","mean"],b=0,lty=2) # dashed line shows overall slope across groups

# Interesting findings:
# The Family Hyaenidae shows a particularly strong relationship between
# mass and premolar length; hyenas are known to crush bones with their teeth (one of
# the functions of the premolar), so it seems understandable that they would have
# an important relationship between a species' mass and the premolar length.
# Ursidae seem to have a relatively weaker relationship between mass and premolar
# length, perhaps because they are omnivorous and don't require large premolars.

