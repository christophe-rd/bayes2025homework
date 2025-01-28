# Homework 1 Bayesian Class
# CRD 15 January 2025

# Goal: complete the workflow for one given dataset. Here I chose penguins

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)
library(palmerpenguins)

fit1 <- FALSE
fit2 <- FALSE

setwd("/Users/christophe_rouleau-desrochers/github/bayes2025homework/homeworksubmitted2/christophe")
# rename dataset
d<-penguins
# === === === === === === === === === === === === === === === === 
#### Step 1. Come up with a model ####
# === === === === === === === === === === === === === === === === 
# Do flippers vary as a function bodymass?
head(d)
str(d)
length(unique(d$species))
adelie <- subset(d, species == "Adelie") # mean flipper length: 190
chinstrap <- subset(d, species == "Chinstrap") # mean flipper length : 195
gentoo <- subset(d, species == "Gentoo") # mean flipper length: 217

# === === === === === === === === === === === === === === === === 
#### Step 2. Simulate data ####
# === === === === === === === === === === === === === === === === 

# Set parameters
nspp <- 3 # nb of species 
mu_flipper <- 100
sigma_flipper <- 10
mu_shift<- 0.1
sigma_shift <- 0.05
species_flipper <- rnorm(nspp, mu_flipper, sigma_flipper)
species_trend <- rnorm(nspp, mu_shift, sigma_shift)
# overall error 
sigma_y <- 50

# Keep the parameters together to compare to model output
paramsgiven <- c(mu_flipper, mu_shift, sigma_shift, sigma_flipper, sigma_y)

# Create the data
mass_min <- 200 # setting it to centigrams

# N per species
n_data_per_species <- round(runif(nspp, 5, 12))
species <- rep(1:nspp, n_data_per_species)
N <- length(species)
mass <- rep(NA, N)

# Simulate mass
### not sure what to do here
mass <- runif(species, min= 200, max=700)
ypred <- length(N)

for (n in 1:N){
  s <- species[n]
  ypred[n] <- species_flipper[s] + species_trend[s]*mass[n]
}

y <- rnorm(N, ypred, sigma_y)

# Plot the data
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(mass), range(y), type="n", xlab="mass", ylab="flipper length",
     bty="l", main="Test data")

for (sp in 1:nspp){
  lines(mass[species==sp], ypred[species==sp], col="darkblue")
}

# Store simulated data in df
df <- data.frame(mass=mass, flipper=ypred, species = species)

# Visualize priors
# Simulate the prior for the intercept
intercept_prior <- rnorm(1000, mean = 100, sd = 20)

# Visualize
pdf("interceptprior_check.pdf", width = 6, height = 8)
hist(intercept_prior, breaks = 30, col = "lightblue", main = "Prior for Intercept",
     xlab = "Flipper Length (Intercept)", probability = TRUE)
abline(v = mean(df$flipper), col = "red", lwd = 2, lty = 2) # Add a line for data mean
dev.off()
# Simulate the prior for the slope
slope_prior <- rnorm(1000, mean = 0.05, sd = 0.02)


# Generate prior predictive values
prior_flipper <- intercept_prior + mean(df$mass) * slope_prior

# Visualize prior predictive vs. actual data
pdf("priorpredictivecheck.pdf", width = 6, height = 8)
hist(prior_flipper, breaks = 30, col = "coral", main = "Prior Predictive Check",
     xlab = "flipper length", probability = TRUE)
abline(v = range(df$flipper), col = "blue", lwd = 2, lty = 2) # Add data range

# Simulate prior for residual standard deviation
residual_prior <- rnorm(1000, mean = 5, sd = 2)  # Example: prior_aux = normal(5, 2)

# Fit the model
fit1 <- stan_lmer(
  flipper ~ mass + (1+mass | species), 
  data = df,                   
  prior_intercept = normal(100, 10),  
  prior = normal(0.05, 0.02),       
  chains = 4,                       
  iter = 2000,                      
  cores = 4                         
)

print(model, digits=2)

# get stan output
sumer <- summary(model)
muparams <- sumer[grep("b", rownames(sumer)), c("mean", "10%", "50%", "90%")]
sigmaparams <- sumer[grep("sigma", rownames(sumer)), c("mean", "10%", "50%", "90%")]

# plot by species levels
spslopes <- sumer[grep("b\\[", rownames(sumer)), "mean"]
spslopes <- sumer[grep("b\\[mass", rownames(sumer)), "mean"]
plot(spslopes~species_trend, xlab="Given species-level slopes", ylab="Modeled species-level slopes", col="darkblue")

### seems to return my slope pretty well, but im not sure

# === === === === === === === === === === === === === === === === 
#### Step 3. Set your priors ####
# === === === === === === === === === === === === === === === === 
# done above

# === === === === === === === === === === === === === === === === 
#### Step 4. Run model on empirical data ####
# === === === === === === === === === === === === === === === === 
head(d)
# Fit the model
fit2 <- stan_lmer(
  flipper_length_mm ~ body_mass_g + (1+body_mass_g | species), 
  data = d,                   
  prior_intercept = normal(100, 10),  
  prior = normal(0.05, 0.02),       
  chains = 4,                       
  iter = 2000,                      
  cores = 4                         
)
print(model)
summary(model)
launch_shinystan(model)
# === === === === === === === === === === === === === === === === 
#### Step 5. Perform retrodictive checks using the model fit to your empiral data ####
# === === === === === === === === === === === === === === === === 
# Posterior predictive checks

yrep <- posterior_predict(model) 

# Predicted vs. Observed
# Calculate the mean predicted values across posterior samples
predicted_mean <- rowMeans(predicted)

# Create a data frame for plotting
predicted_vs_observed <- data.frame(
  Observed = d$flipper_length_mm,   # Observed values
  Predicted = predicted_mean       # Mean predicted values
)
summary(model)

# Extract predictions
posterior_pred <- as.data.frame(posterior_predict(model))
rmvrows<-d[1:342,] # removing 2 rows because not same nb of rows
rmvrows$predicted <- apply(posterior_pred, 2, mean)  # Mean predicted value for each observation


# Observed vs Predicted Flipper by Species
postpredict<-ggplot(rmvrows, aes(x = flipper_length_mm, y = predicted, color = factor(species))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "observed vs predicted",
       x = "Observed flipper length",
       y = "Predicted flipper length",
       color = "Species") +
  theme_minimal()
ggsave("posteriorcheck.pdf", postpredict)

### This doesn't make sens as for 2 species, the predicted value decreases as the observed value increases. 
### I should have spent more time working on my priors. Should have started working on this earlier!
