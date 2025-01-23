#install.packages("rstanarm")
## Jan 2025 Wenhuan Xu

library(rstanarm)

#set.seed(123)

# set sample size
n <- 100
a <- 1 # intercept

#set a continuous variable
x1 <-runif(n, min=0, max=100)
b1 <- 6
hist(x1)

#set a binary covariate 
x2 <- sample(c(0,1), n, replace=TRUE)
b2 <- 4
b3 <- 2 # interaction effect size

#set error term
#error <- rnorm(n, mean=0, sd=3)
#error <- rnorm(n, mean=0, sd=18)
error <- rnorm(n, mean=0, sd=9)


#create similuated data
y <- a+b1*x1 +b2*x2 +b3*x1*x2 +error
hist(y)
data <- data.frame(x1=x1, x2=x2, y=y)
head(data)

# visualize similated data
plot(x1,y, main="Fish weight and body size", xlab="x", ylab="y", pch=16, col="blue")
abline(a=1, b=6, col="red", lwd=2)
plot(x2,y, main="Fish weight and treatment", xlab="x", ylab="y", pch=16, col="blue")

m0<-lm(y~x1*x2) ## bayesian model
summary(m0)
coef(m0)

m1<-stan_glm(y~x1*x2) ## bayesian model
summary(m1)
coef(m1)
#(Intercept)          x1          x2       x1:x2 
#  0.5005876   6.0065799   4.0990015   1.9928442
# Answer:generally good, but intercept is deviated


## change sample size to 20%
set.seed(123)
data20_indice <- sample(1:n, size =20, replace = FALSE)
data20 <- data[data20_indice,]
head(data20)
m_20 <- stan_glm(y~x1*x2, data =data20)
coef(m_20)
#(Intercept)          x1          x2       x1:x2 
#-0.2116596   6.0244183   3.8180087   1.9966082
# Answer: the intercept is deviated a lot, effect size of x2,interation sometimes are not good.

n<-100
### sampling from 1-n and plot it
n_samples <- 5:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = 4)
colnames(bayes_coefficients) <- names(coef(m0))

for (i in n_samples) {
  sample_data <- data[1:i, ]
  m2_sample <- stan_glm(y ~ x1 * x2, data = sample_data)
  bayes_coefficients[i, ] <- coef(m2_sample)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_sd18_100.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(2, 2)) # Arrange plots

for (j in 1:4) {
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2) # True parameter value
}
dev.off()
par(mfrow = c(1, 1)) # Reset plotting layout



### improved sampling from 1-n with 10 draws each time
n_samples <- 5:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = 4)
colnames(bayes_coefficients) <- names(coef(m0))

for (i in n_samples) {
  sample_data <- data[1:i, ]
  draws <- replicate(10, {
    sample_indices <- sample(1:n, size = i, replace = FALSE) # Randomly select 'i' rows
    sample_data <- data[sample_indices, ]
    m2_sample  <- stan_glm(y ~ x1 * x2, data = sample_data)
    coef(m2_sample)
    
  })
  bayes_coefficients[i, ] <- rowMeans(draws)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_sd18_100_10rep.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(2, 2)) # Arrange plots
for (j in 1:4) {
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2) # True parameter value
}
dev.off()
par(mfrow = c(1, 1)) # Reset plotting layout.






## Step 3
#change error term
n <- 100
error <- rnorm(n, mean=0, sd=9)

#create similuated data
y <- a+b1*x1 +b2*x2 +b3*x1*x2 +error
data <- data.frame(x1=x1, x2=x2, y=y)
head(data)

# visualize similated data
plot(x1,y, main="Fish weight and body size", xlab="x", ylab="y", pch=16, col="blue")
abline(a=1, b=6, col="red", lwd=2)
plot(x2,y, main="Fish weight and treatment", xlab="x", ylab="y", pch=16, col="blue")

m1<-stan_glm(y~x1*x2) ## bayesian model
summary(m1)
coef(m1)
#(Intercept)          x1          x2       x1:x2 
#3.591158    5.977960    1.317760    2.020746


## change sample size to 20%
set.seed(123)
data20_indice <- sample(1:n, size =20, replace = FALSE)
data20 <- data[data20_indice,]
head(data20)
m_20 <- stan_glm(y~x1*x2, data =data20)
coef(m_20)
#(Intercept)          x1          x2       x1:x2 
#18.346338    4.122492    3.757669   14.087005 

n<-20
### sampling from 1-n and plot it
n_samples <- 2:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = 4)
colnames(bayes_coefficients) <- names(coef(m0))

for (i in n_samples) {
  sample_data <- data[1:i, ]
  m2_sample <- stan_glm(y ~ x1 * x2, data = sample_data)
  bayes_coefficients[i, ] <- coef(m2_sample)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_sd9.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(2, 2)) # Arrange plots
for (j in 1:4) {
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2) # True parameter value
}
dev.off()
par(mfrow = c(1, 1)) # Reset plotting layout



### improved sampling from 1-n with 10 draws each time
n_samples <- 2:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = 4)
colnames(bayes_coefficients) <- names(coef(m0))

for (i in n_samples) {
  sample_data <- data[1:i, ]
  draws <- replicate(10, {
    sample_indices <- sample(1:n, size = i, replace = FALSE) # Randomly select 'i' rows
    sample_data <- data[sample_indices, ]
    m_sample <- stan_glm(y ~ x1 * x2, data = sample_data)
    coef(m_sample)
  })
  bayes_coefficients[i, ] <- rowMeans(draws)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_sd9_rep10.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(2, 2)) # Arrange plots
for (j in 1:4) {
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2) # True parameter value
}
dev.off()
par(mfrow = c(1, 1)) # Reset plotting layout.





## Chanlenge C

# Fit Bayesian model to the dataset
m1 <- stan_glm(y ~ x1*x2, data = data)

# Extract 50% uncertainty intervals for all parameters
posterior_summary <- posterior_interval(m1, prob = 0.50)

# Extract posterior medians
posterior_medians <- coef(m1)

# Check parameter names
summary_param_names <- rownames(posterior_summary)
median_param_names <- names(posterior_medians)

# Find matching parameters
common_params <- intersect(summary_param_names, median_param_names)

# Subset to only matching parameters
posterior_summary <- posterior_summary[common_params, , drop = FALSE]
posterior_medians <- posterior_medians[common_params]

# Combine results into a data frame
results <- data.frame(
  Parameter = common_params,
  Median = posterior_medians,
  Lower = posterior_summary[, 1],
  Upper = posterior_summary[, 2]
)

# Print results to verify
print(results)

# Define sample sizes for the loop
sample_sizes <- seq(10, 90, by = 10)

# Initialize a list to store results for each sample size
results_list <- list()

# Loop through each sample size
for (i in seq_along(sample_sizes)) {
  sample_size <- sample_sizes[i]
  
  # Subset data
  sample_data <- data[1:sample_size, ]
  
  # Fit Bayesian model
  model <- stan_glm(y ~ x1*x2, data = sample_data)
  
  # Extract 50% uncertainty intervals
  posterior_summary <- posterior_interval(model, prob = 0.50)
  
  # Extract posterior medians
  posterior_medians <- coef(model)
  
  # Find matching parameters
  summary_param_names <- rownames(posterior_summary)
  median_param_names <- names(posterior_medians)
  common_params <- intersect(summary_param_names, median_param_names)
  
  # Subset to only matching parameters
  posterior_summary <- posterior_summary[common_params, , drop = FALSE]
  posterior_medians <- posterior_medians[common_params]
  
  # Combine results into a data frame
  results <- data.frame(
    Parameter = common_params,
    Median = posterior_medians,
    Lower = posterior_summary[, 1],
    Upper = posterior_summary[, 2]
  )
  
  # Store results
  results_list[[i]] <- list(
    sample_size = sample_size,
    results = results
  )
}

# Create a 3x3 grid of plots
tiff(filename = "uncertainty_intervals_sd9_9plots.tiff", width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(5, 8, 4, 2)) # 3x3 grid and adjusted margins

# Plot results for each sample size
for (i in seq_along(results_list)) {
  sample_size <- results_list[[i]]$sample_size
  results <- results_list[[i]]$results
  
  # Plot posterior estimates with 50% uncertainty intervals
  plot(
    x = results$Median, 
    y = seq_along(results$Parameter), 
    xlim = range(c(results$Lower, results$Upper)), 
    ylim = c(0.5, length(results$Parameter) + 0.5),
    pch = 16, col = "blue", xlab = "Estimate", ylab = "",
    main = paste("Sample Size:", sample_size), axes = FALSE
  )
  
  # Add intervals
  segments(
    x0 = results$Lower, x1 = results$Upper,
    y0 = seq_along(results$Parameter), y1 = seq_along(results$Parameter),
    col = "blue", lwd = 2
  )
  
  # Customize axes
  axis(1) # X-axis for estimates
  axis(2, at = seq_along(results$Parameter), labels = results$Parameter, las = 2) # Y-axis for parameters
  
  # Add a reference line at 0
  abline(v = 0, col = "red", lty = 2)
}

dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))






## ten draws
# Define sample sizes for the loop
sample_sizes <- seq(10, 90, by = 10)

# Initialize a list to store results for each sample size
results_list <- list()

# Loop through each sample size
for (i in seq_along(sample_sizes)) {
  sample_size <- sample_sizes[i]
  
  # Store results for 10 repetitions
  draw_results <- replicate(10, {
    # Randomly sample rows
    sample_indices <- sample(1:n, size = sample_size, replace = FALSE)
    sample_data <- data[sample_indices, ]
    
    # Fit Bayesian model
    model <- stan_glm(y ~ x1 + x2 + x1:x2, data = sample_data)
    
    # Extract posterior medians and 50% uncertainty intervals
    posterior_summary <- posterior_interval(model, prob = 0.50)
    posterior_medians <- coef(model)
    
    # Find matching parameters
    summary_param_names <- rownames(posterior_summary)
    median_param_names <- names(posterior_medians)
    common_params <- intersect(summary_param_names, median_param_names)
    
    # Subset to only matching parameters
    posterior_summary <- posterior_summary[common_params, , drop = FALSE]
    posterior_medians <- posterior_medians[common_params]
    
    # Combine results into a data frame
    data.frame(
      Parameter = common_params,
      Median = posterior_medians,
      Lower = posterior_summary[, 1],
      Upper = posterior_summary[, 2]
    )
  }, simplify = FALSE)
  
  # Average results across 10 repetitions
  combined_results <- do.call(rbind, draw_results)
  averaged_results <- aggregate(. ~ Parameter, data = combined_results, FUN = mean)
  
  # Store averaged results
  results_list[[i]] <- list(
    sample_size = sample_size,
    results = averaged_results
  )
}

# Create a 3x3 grid of plots
tiff(filename = "uncertainty_intervals_10reps_9plots.tiff", width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(5, 8, 4, 2)) # 3x3 grid and adjusted margins

# Plot results for each sample size
for (i in seq_along(results_list)) {
  sample_size <- results_list[[i]]$sample_size
  results <- results_list[[i]]$results
  
  # Plot posterior estimates with 50% uncertainty intervals
  plot(
    x = results$Median, 
    y = seq_along(results$Parameter), 
    xlim = range(c(results$Lower, results$Upper)), 
    ylim = c(0.5, length(results$Parameter) + 0.5),
    pch = 16, col = "blue", xlab = "Estimate", ylab = "",
    main = paste("Sample Size:", sample_size), axes = FALSE
  )
  
  # Add intervals
  segments(
    x0 = results$Lower, x1 = results$Upper,
    y0 = seq_along(results$Parameter), y1 = seq_along(results$Parameter),
    col = "blue", lwd = 2
  )
  
  # Customize axes
  axis(1) # X-axis for estimates
  axis(2, at = seq_along(results$Parameter), labels = results$Parameter, las = 2) # Y-axis for parameters
  
  # Add a reference line at 0
  abline(v = 0, col = "red", lty = 2)
}

dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))
