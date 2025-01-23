library(rstanarm)

set.seed(123)

# Set sample size
n <- 100
a <- 1 # Intercept

# Create a continuous variable
x1 <- runif(n, min = 0, max = 100)
b1 <- 6

# Create treatment variables
treatment1 <- factor(sample(c("Control", "Treatment1"), n, replace = TRUE))
treatment2 <- factor(sample(c("Control", "Treatment2"), n, replace = TRUE))

# Manually create dummy variables
treatment1_Treatment <- ifelse(treatment1 == "Treatment1", 1, 0)
treatment2_Treatment <- ifelse(treatment2 == "Treatment2", 1, 0)

# Set coefficients for the dummy variables and interactions
b_treatment1 <- 4
b_treatment2 <- 3
b_interaction_treatment1 <- 2
b_interaction_treatment2 <- -1
b_interaction_treatment1_treatment2 <- 1.5
b_three_way_interaction <- -0.5

# Set error term
error <- rnorm(n, mean = 0, sd = 18)

# Simulate response variable with all interactions
y <- a + b1 * x1 +
  b_treatment1 * treatment1_Treatment +
  b_treatment2 * treatment2_Treatment +
  b_interaction_treatment1 * x1 * treatment1_Treatment +
  b_interaction_treatment2 * x1 * treatment2_Treatment +
  b_interaction_treatment1_treatment2 * treatment1_Treatment * treatment2_Treatment +
  b_three_way_interaction * x1 * treatment1_Treatment * treatment2_Treatment +
  error

# Create data frame with dummy variables
data <- data.frame(x1 = x1, 
                   treatment1_Treatment = treatment1_Treatment, 
                   treatment2_Treatment = treatment2_Treatment, 
                   y = y)
head(data)

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
  model <- stan_glm(y ~ x1 + treatment1_Treatment + treatment2_Treatment + 
                      x1:treatment1_Treatment + x1:treatment2_Treatment + 
                      treatment1_Treatment:treatment2_Treatment + 
                      x1:treatment1_Treatment:treatment2_Treatment, 
                    data = sample_data)
  
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
tiff(filename = "uncertainty_intervals_full_factorial_9plots.tiff", width = 14, height = 10, units = "in", res = 300)
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
