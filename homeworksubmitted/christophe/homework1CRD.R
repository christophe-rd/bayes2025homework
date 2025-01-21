# Homework 1 Bayesian Class
# CRD 15 January 2025

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)

# === === === === === === === === === === === === === === === === 
#### Step 1. Simulate data ####
# === === === === === === === === === === === === === === === === 

### My linear regression will  analyse how growth increment (mm) is affected by the number of days in the growing season
### Predictor variable : number of days
### Response variable : growth in increment in ring width (mm)
n <- 100
a <- 10 
b <- 0.8
sigma <- 15 

# Simulate x
x <- runif(n, min=100, max=150)
# Generate error with a mean of 0
error<-rnorm(n, 0, sigma)
# Z-score standardization of x (number of days)
x_z <- (x - mean(x)) / (2*sd(x))
# Z-score standardization of error
error_z <- (error - mean(error)) / (2*sd(error))
# Simulate y
ypred <- a+b*x_z+error_z
# Plot y ~ x
plot(ypred~x)
dev.off()

# Fit the model using stan glm
fit_1 <- stan_glm(ypred ~ x_z)
# Display the results
print(fit_1, digits=2) 
# Answer: the model returned 0.7 for the slope and 9.76 for the intercept which is very close to the parameters I set (0.8 and 10 respectively)

# Plot the data and fitted regression line
plot(x_z, ypred, main="Data and fitted regression line") 
a_hat <- coef(fit_1)[1] 
b_hat <- coef(fit_1)[2] 
abline(a_hat, b_hat)
x_bar <- mean(x) 
text(x_bar, a_hat + b_hat*x_bar, paste("y =", round(a_hat, 2), "+", round(b_hat, 2), "* x"), adj=c(3,-2))
dev.off()

##### Add binary covariate #####
# Simulate treatment (binary covariate)
treatment <- rbinom(n, 1, prob = 0.5) 
# Setting effect sizes:
c <- 0.1      # Effect size of treatment which has lower effect size than my predictor
d <- b/2        # Effect size of interaction : half of b
# Simulate Y (growth increment) with treatment
ybin <- a + b*x_z + c*treatment + d*(x_z * treatment) + error_z
# Plot
plot(x_z, ybin, col=treatment + 1, pch=20,
     main="", xlab="Days in Growing Season", ylab="Growth Increment")
legend("topright", legend=c("Control", "Treatment"), col=c(1, 2), pch=16)
# Df
simbin <- data.frame(x_z, treatment, ybin)

fit_2 <- stan_glm(ybin ~ x_z*treatment, data = simbin, cores=4)
print(fit_2, digits=2) 
# Answer: I return the parameter value not as well as without the treatment (number of days:0.73 compared to 0.8, treatment: 0.08 compared to 0.4 and x_z:treatment: 0.42 compared to 0.2)

# === === === === === === === === === === === === === === === === 
#### Step 2 ####
# === === === === === === === === === === === === === === === === 
##### Using 20% of my data ####
# Modify parameters with 20 suffix
n20 <- n / 5
# Simulate x
x20 <- runif(n20, min=100, max=150)
# Generate error with a mean of 0
error20 <- rnorm(n20, 0, sigma)
# Treatment
treatment20 <- rbinom(n20, 1, prob=0.5)
# Z-score standardization of x (number of days)
x_z20 <- (x20 - mean(x20)) / (2*sd(x20))
# Z-score standardization of error
error_z20 <- (error20 - mean(error20)) / (2*sd(error20))
# Simulate Y (growth increment) with treatment and standardized x
ybin20 <- a + b * x_z20 + c * treatment20 + d * (x_z20 * treatment20) + error_z20
# Plot simulated data with n of 20
plot(x_z20, ybin20, col=treatment20 + 1, pch=20,
     main="", xlab="Days in Growing Season", ylab="Growth Increment")
legend("topright", legend=c("Control", "Treatment"), col=c(1, 2), pch=16)
# Create data frame
simbin20 <- data.frame(x_z20, treatment20, ybin20)
# Fit the model
fit_3 <- stan_glm(ybin20 ~ x_z20 * treatment20, data = simbin20, cores=4)
# Print results
print(fit_3, digits=2)
# Answer: I don't return the parameters very well (number of days: 0.71 compared to 0.8, treatment: -0.07 compared to 0.4 and interactioln: 0.34 compared to 0.2). I don't understand why the effect size of the treatment is now negative. 

# Create a data frame for plotting data from 1.
samplesize <- c(20, 100)
x <- c(coef(fit_3)[2], coef(fit_2)[2])
treatments <- c(coef(fit_3)[3], coef(fit_2)[3])
interaction <- c(coef(fit_3)[4], coef(fit_2)[4])
df4plot <- data.frame(samplesize, x = x, treatments = treatments, interaction = interaction)
# Plot of sampling from my data and parameter comparison BEFORE improving sampling
ggplot(df4plot) +
  geom_point(aes(x=samplesize, y=x), color="red") +
  geom_point(aes(x=samplesize, y=treatments), color="green") +
  geom_point(aes(x=samplesize, y=interaction), color="blue") +
  scale_x_continuous(breaks = c(20, 100), limits = c(0, 120))+   labs(
    x = "Sample Size",
    y = "Values",
    title = "Estimated Parameters by Sample Size",
    subtitle = "Red: number of days, Green: Treatment, Blue: Nb of days X Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),  
    legend.position = "Left"  
  )

##### Improve sampling ##### 
# store results
results_df <- data.frame(iteration = integer(), estimate = numeric(), stringsAsFactors = FALSE)
# Loop for n = 100 and n=20
results_df <- data.frame(iteration = integer(), n = integer(), intercept = numeric(),
                         coef_x = numeric(), coef_treatment = numeric(), coef_interaction = numeric(), 
                         stringsAsFactors = FALSE)
n <- 100
for (i in 1:10) {  # repeat 10 times for 10 draws
  # Simulate data
  x <- runif(n, min = 100, max = 150)  
  error <- rnorm(n, 0, sigma)          
  treatment <- rbinom(n, 1, prob = 0.5)
  # z score with 2 sd
  x_z <- (x - mean(x)) / (2 * sd(x))
  error_z <- (error - mean(error)) / (2 * sd(error))  
  ybin <- a + b*x_z + c*treatment + d*(x_z*treatment) + error_z  
  
  # Create the simulated dataset
  simbin <- data.frame(ybin = ybin, x = x_z, treatment = treatment)  
  # fit the model
  fitloop <- stan_glm(ybin ~ x * treatment, data = simbin)  
  # Extract the coefs
  coefs <- coef(fitloop)
  # Store results
  results_df <- rbind(results_df, data.frame(
    iteration = i, 
    n = 100,
    intercept = coefs["(Intercept)"], 
    coef_x = coefs["x"], 
    coef_treatment = coefs["treatment"], 
    coef_interaction = coefs["x:treatment"]
  ))
}

# Loop for n = 20
n <- 20
for (i in 1:10) {  # repeat 10 times for 10 draws
  # Simulate data
  x <- runif(n, min = 100, max = 150)  
  error <- rnorm(n, 0, sigma)          
  treatment <- rbinom(n, 1, prob = 0.5)
  # z score with 2 sd
  x_z <- (x - mean(x)) / (2 * sd(x))
  error_z <- (error - mean(error)) / (2 * sd(error))  
  ybin <- a + b*x_z + c*treatment + d*(x_z*treatment) + error_z  
  simbin <- data.frame(ybin = ybin, x = x_z, treatment = treatment)  
  # fit the model
  fitloop <- stan_glm(ybin ~ x * treatment, data = simbin)
  # Extract coefs
  coefs <- coef(fitloop)
  # Store the results for n = 20
  results_df <- rbind(results_df, data.frame(
    iteration = i, 
    n = 20,
    intercept = coefs["(Intercept)"], 
    coef_x = coefs["x"], 
    coef_treatment = coefs["treatment"], 
    coef_interaction = coefs["x:treatment"]
  ))
}

# calculate mean and SE for plot
unique_n <- unique(results_df$n)
mean_x <- numeric(length(unique_n))
se_x <- numeric(length(unique_n))

# get the mean and sd of the iterations
for (i in 1:length(unique_n)) {
  xforn <- results_df$coef_x[results_df$n == unique_n[i]]
  mean_x[i] <- mean(xforn)
  se_x[i] <- sd(xforn) / sqrt(length(xforn))
}

# Create df
plot_data <- data.frame(n = unique_n, mean_x = mean_x, se_x = se_x)
# Plot
ggplot(plot_data, aes(x = n, y = mean_x)) +
  geom_point(size = 3, color = "blue") + 
  geom_errorbar(aes(ymin = mean_x - se_x, ymax = mean_x + se_x),
                width = 2) +
  labs(x = "n", y = "Mean Intercept ± SE") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  ggtitle("Mean nb of days with SE with 10 iterations") +
  scale_x_continuous(breaks = c(20, 100), limits = c(0, 120))+   labs(
    x = "Sample Size",
    y = "Parameter estimates",
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),  
    legend.position = "Left" )

### didn't have time to add the other parameters
# Answer: for this parameter only (days in the GS), there is far more variability for the estimated parameters of n20 compared to n100.

# === === === === === === === === === === === === === === === === 
#### Step 3. Change error term by 50% ####
# === === === === === === === === === === === === === === === === 
# store results
results_df <- data.frame(iteration = integer(), estimate = numeric(), stringsAsFactors = FALSE)
# Loop for n = 100 and n=20
results_df <- data.frame(iteration = integer(), n = integer(), intercept = numeric(),
                         coef_x = numeric(), coef_treatment = numeric(), coef_interaction = numeric(), 
                         stringsAsFactors = FALSE)
n <- 100
for (i in 1:10) {  # repeat 10 times for 10 draws
  # Simulate data
  x <- runif(n, min = 100, max = 150)  
  error <- rnorm(n, 0, sigma)*0.5          
  treatment <- rbinom(n, 1, prob = 0.5)
  # z score with 2 sd
  x_z <- (x - mean(x)) / (2 * sd(x))
  error_z <- (error - mean(error)) / (2 * sd(error))  
  ybin <- a + b*x_z + c*treatment + d*(x_z*treatment) + error_z  
  
  # Create the simulated dataset
  simbin <- data.frame(ybin = ybin, x = x_z, treatment = treatment)  
  # fit the model
  fitloop <- stan_glm(ybin ~ x * treatment, data = simbin)  
  # Extract the coefs
  coefs <- coef(fitloop)
  # Store results
  results_df <- rbind(results_df, data.frame(
    iteration = i, 
    n = 100,
    intercept = coefs["(Intercept)"], 
    coef_x = coefs["x"], 
    coef_treatment = coefs["treatment"], 
    coef_interaction = coefs["x:treatment"]
  ))
}

# Loop for n = 20
n <- 20
for (i in 1:10) {  # repeat 10 times for 10 draws
  # Simulate data
  x <- runif(n, min = 100, max = 150)  
  error <- rnorm(n, 0, sigma)*0.5          
  treatment <- rbinom(n, 1, prob = 0.5)
  # z score with 2 sd
  x_z <- (x - mean(x)) / (2 * sd(x))
  error_z <- (error - mean(error)) / (2 * sd(error))  
  ybin <- a + b*x_z + c*treatment + d*(x_z*treatment) + error_z  
  simbin <- data.frame(ybin = ybin, x = x_z, treatment = treatment)  
  # fit the model
  fitloop <- stan_glm(ybin ~ x * treatment, data = simbin)
  # Extract coefs
  coefs <- coef(fitloop)
  # Store the results for n = 20
  results_df <- rbind(results_df, data.frame(
    iteration = i, 
    n = 20,
    intercept = coefs["(Intercept)"], 
    coef_x = coefs["x"], 
    coef_treatment = coefs["treatment"], 
    coef_interaction = coefs["x:treatment"]
  ))
}

# calculate mean and SE for plot
unique_n <- unique(results_df$n)
mean_x <- numeric(length(unique_n))
se_x <- numeric(length(unique_n))

# get the mean and sd of the iterations
for (i in 1:length(unique_n)) {
  xforn <- results_df$coef_x[results_df$n == unique_n[i]]
  mean_x[i] <- mean(xforn)
  se_x[i] <- sd(xforn) / sqrt(length(xforn))
}

# Create df
plot_data <- data.frame(n = unique_n, mean_x = mean_x, se_x = se_x)
# Plot
gg<-ggplot(plot_data, aes(x = n, y = mean_x)) +
  geom_point(size = 3, color = "blue") + 
  geom_errorbar(aes(ymin = mean_x - se_x, ymax = mean_x + se_x),
                width = 2) +
  labs(x = "n", y = "Mean Intercept ± SE") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  ggtitle("Mean nb of days with SE with 10 iterations") +
  scale_x_continuous(breaks = c(20, 100), limits = c(0, 120))+   labs(
    x = "Sample Size",
    y = "Parameter estimates",
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),  
    legend.position = "Left" )

# Answer: for the parameter x only, it seems like reducing the error by 50% lowered the SE a little. 

### didn't have time to compare the effect of a reduced error accross the different parameter values and sample sizes... 
