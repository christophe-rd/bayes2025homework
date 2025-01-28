#library(truncnorm)
library(rstan)



# Step1 is the stan file named partial_pooling.stan



# Step2 simulate test data
# Set seed for reproducibility
set.seed(42)

# Number of species
n_species <- 10

# Generate sample sizes for each species (varying)
sample_sizes <- sample(20:100, n_species, replace = TRUE)

# Global-level parameters
mu_intercept <- 10  # Mean intercept
mu_slope <- 0.8       # Mean slope
sigma_intercept <- 5  # Standard deviation of intercepts
sigma_slope <- 0.2    # Standard deviation of slopes
sigma_residual <- 3   # Residual standard deviation

# Generate species-specific intercepts and slopes
species_intercepts <- rnorm(n_species, mean = mu_intercept, sd = sigma_intercept)
species_slopes <- rnorm(n_species, mean = mu_slope, sd = sigma_slope)

# Initialize empty data frame
stan_data <- data.frame(
  species = integer(),
  ind = numeric(),
  y = numeric()
)

# Generate data for each species
for (i in 1:n_species) {
  # Number of samples for this species
  n_samples <- sample_sizes[i]
  
  # Predictor variable (e.g., body length)
  ind <- rnorm(n_samples, mean = 30, sd = 5)
  
  # Response variable (e.g., body size)
  y <- species_intercepts[i] + species_slopes[i] * ind + rnorm(n_samples, mean = 0, sd = sigma_residual)
  
  # Create species-specific data
  species_data <- data.frame(
    species = rep(i, n_samples),  # Species as integer IDs
    ind = ind,
    y = y
  )
  
  # Combine into the main data frame
  stan_data <- rbind(stan_data, species_data)
}

# Check the first few rows of the data
head(stan_data)
hist(y)
# Save the data for Stan
#write.csv(stan_data, "stan_model_data.csv", row.names = FALSE)

# Load ggplot2 library
library(ggplot2)

# Convert species to a factor for better visualization
stan_data$species <- as.factor(stan_data$species)

# Plot with ggplot
simulated_plot<-ggplot(stan_data, aes(x = ind, y = y, color = species)) +
  geom_point(alpha = 0.6) +  # Scatter points with transparency
  geom_smooth(method = "lm", se = FALSE, aes(group = species)) +  # Linear trend lines for each species
  labs(
    title = "Scatter Plot with Trend Lines for Each Species",
    x = "Ind (Body Length)",
    y = "Y (Body Size)",
    color = "Species"
  ) +
  theme_minimal() +  # Minimal theme for better aesthetics
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )
tiff(filename = "simulated_data.tiff", width = 10, height = 7, units = "in", res = 300)
simulated_plot
dev.off()

# data list
stan_data_list <- list(
  N = nrow(stan_data),                       # Total number of observations
  Nspp = length(unique(stan_data$species)),  # Number of species (groups)
  species = as.integer(stan_data$species),   # Species IDs as integers
  ind = stan_data$ind,                     # Predictor variable (year or body length)
  y = stan_data$y                            # Response variable (body size)
)

# Verify the structure 
str(stan_data_list)

# Compile the Stan model
stan_file_path <- "D:/Post-Doc information/Hyierecial Bayesian workshop_Ellize/Assignment2/carnivores/partial_pooling.stan"
stan_model <- stan_model(file = stan_file_path)

# Fit the model
fit <- sampling(
  object = stan_model,
  data = stan_data_list,
  iter = 2000,    # Number of iterations
  chains = 4,     # Number of chains
  seed = 42       # Seed for reproducibility
)

# Print the summary of the fit
print(fit)


# stan output
sumer<- summary(fit)$summary
muparams <- sumer[grep("mu",rownames(sumer)), c("mean", "25%", "50%", "75%")]
sigmaparams <- sumer[grep("sigma", rownames(sumer)), c("mean", "25%", "50%", "75%")]
paramsgiven <- c(mu_intercept, mu_slope, sigma_intercept, sigma_slope, sigma_residual)
paramsgiven
muparams
sigmaparams

spslopes <- sumer[grep("b\\[", rownames(sumer)), "mean"]

plot(spslopes~species_slopes, xlab="Given species-level slopes", ylab="Modeled species-level slopes", col="darkblue")
abline(0,1)

## the result is good 




## Step 4 use empirical data
# Load the data (assuming the file is named 'carnivoreteeth.csv' and is in your working directory)
carnivore_data <- read.csv("carnivoreteeth.csv", header= TRUE)
head(carnivore_data)
colnames(carnivore_data)[1] <- "Species"
data0<-subset(carnivore_data, Lat>0)
data1<- na.omit(subset(data0, sex=="Male", select=c(Family, PM4, CsupL)))
head(data1)

data1$Family <- as.integer(as.factor(data1$Family))

# data list
real_data_list <- list(
  N = nrow(data1),                       # Total number of observations
  Nspp = length(unique(data1$Family)),  # Number of species (groups)
  species = data1$Family,   # Species IDs as integers
  ind = data1$PM4,                     # Predictor variable (year or body length)
  y = data1$CsupL                            # Response variable (body size)
)

# Verify the structure 
str(real_data_list)

# Compile the Stan model
stan_file_path2 <- "D:/Post-Doc information/Hyierecial Bayesian workshop_Ellize/Assignment2/carnivores/partial_pooling.stan"
stan_model2 <- stan_model(file = stan_file_path2)

# Fit the model
fit2 <- sampling(
  object = stan_model2,
  data = real_data_list,
  iter = 2000,    # Number of iterations
  chains = 4,     # Number of chains
  seed = 42       # Seed for reproducibility
)

# Print the summary of the fit
print(fit2)


# stan output
sumer2<- summary(fit2)$summary
muparams2 <- sumer2[grep("mu",rownames(sumer2)), c("mean", "25%", "50%", "75%")]
sigmaparams2 <- sumer2[grep("sigma", rownames(sumer2)), c("mean", "25%", "50%", "75%")]
muparams2
sigmaparams2

spslopes <- sumer[grep("b\\[", rownames(sumer)), "mean"]

plot(spslopes~species_slopes, xlab="Given species-level slopes", ylab="Modeled species-level slopes", col="darkblue")
abline(0,1)






# Step 5: do retrodictive checks
posterior_samples <- extract(fit2)

# Extract the slopes (b[Nspp])
slopes <- posterior_samples$b  # b is a matrix with Nspp rows (species) and posterior draws


# Calculate summary statistics for each species slope
slopes_summary <- apply(slopes, 2, function(x) {
  c(mean = mean(x), 
    lower_95 = quantile(x, 0.025), 
    upper_95 = quantile(x, 0.975))
})

# Ensure slopes_summary is transposed correctly
slopes_summary <- t(slopes_summary)


slopes_df <- data.frame(
  Species = 1:nrow(slopes_summary),  # Species ID
  Mean = slopes_summary[, "mean"],   # Mean of slopes
  Lower_95 = slopes_summary[, "lower_95.2.5%"],  # Lower 95% CI
  Upper_95 = slopes_summary[, "upper_95.97.5%"]   # Upper 95% CI
)

print(slopes_df)

head(posterior_samples$a)

# Initialize the posterior_predictions matrix
posterior_predictions <- matrix(NA, nrow = nrow(slopes), ncol = real_data_list$N)

# Generate posterior predictive samples
for (i in 1:nrow(slopes)) { # Loop over posterior iterations
  for (j in 1:real_data_list$N) { # Loop over data points
    posterior_predictions[i, j] <- 
      posterior_samples$a[i, real_data_list$species[j]] + 
      slopes[i, real_data_list$species[j]] * real_data_list$ind[j]
  }
}


# Calculate the mean and credible intervals for the predictions
pred_mean <- apply(posterior_predictions, 2, mean)
pred_lower <- apply(posterior_predictions, 2, quantile, 0.025)
pred_upper <- apply(posterior_predictions, 2, quantile, 0.975)

# Plot observed vs. predicted
library(ggplot2)
ppc_data <- data.frame(
  Observed = real_data_list$y,
  Predicted = pred_mean,
  Lower_95 = pred_lower,
  Upper_95 = pred_upper
)

ppcfig<-ggplot(ppc_data, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_errorbar(aes(ymin = Lower_95, ymax = Upper_95), width = 0.2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Posterior Predictive Check",
    x = "Observed Response (y)",
    y = "Predicted Response (Posterior Mean)"
  ) +
  theme_minimal()

tiff(filename = "PosteriorPredCheck.tiff", width = 10, height = 7, units = "in", res = 300)
ppcfig
dev.off()

slope_fig<- ggplot(slopes_df, aes(x = Species, y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower_95, ymax = Upper_95), width = 0.2) +
  labs(
    title = "Posterior Distribution of Slopes by Species",
    x = "Species ID",
    y = "Slope (Posterior Mean)"
  ) +
  theme_minimal()
tiff(filename = "PosteriorPredCheck_slope.tiff", width = 10, height = 7, units = "in", res = 300)
slope_fig
dev.off()



library(shinystan)
# Launch the shinystan interface for the Stan model
launch_shinystan(fit2)
