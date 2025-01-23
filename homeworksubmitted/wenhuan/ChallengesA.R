library(rstanarm)

set.seed(123)

# Set sample size
n <- 100
a <- 1 # Intercept

# Create a continuous variable
x1 <- runif(n, min = 0, max = 100)
b1 <- 6

# Create treatment variable with 3 levels
x2 <- factor(sample(c("Control", "Treatment1", "Treatment2"), n, replace = TRUE))
# Manually create dummy variables for Treatment1 and Treatment2 (Control as baseline)
x2_Treatment1 <- ifelse(x2 == "Treatment1", 1, 0)
x2_Treatment2 <- ifelse(x2 == "Treatment2", 1, 0)

# Set coefficients for the dummy variables and interactions
b_Treatment1 <- 4
b_Treatment2 <- 2
b_interaction_Treatment1 <- 2
b_interaction_Treatment2 <- -1

# Set error term
#error <- rnorm(n, mean = 0, sd = 18)
error <- rnorm(n, mean = 0, sd = 9)

# Simulate response variable
y <- a + b1 * x1 +
  b_Treatment1 * x2_Treatment1 +
  b_Treatment2 * x2_Treatment2 +
  b_interaction_Treatment1 * x1 * x2_Treatment1 +
  b_interaction_Treatment2 * x1 * x2_Treatment2 +
  error

# Create data frame with dummy variables
data <- data.frame(x1 = x1, x2_Treatment1 = x2_Treatment1, x2_Treatment2 = x2_Treatment2, y = y)
head(data)

# Fit a Bayesian model with dummy variables
m1 <- stan_glm(y ~ x1 + x2_Treatment1 + x2_Treatment2 + 
                 x1:x2_Treatment1 + x1:x2_Treatment2, data = data)
summary(m1)
coef(m1)
#(Intercept)    x1     x2_Treatment1    x2_Treatment2  x1:x2_Treatment1 x1:x2_Treatment2 
#-9.613722    6.179438   4.375515         6.200419      1.953777        -1.112037 
#
#


# Repeat for a smaller sample size (20%)
set.seed(123)
data20_indice <- sample(1:n, size = 20, replace = FALSE)
data20 <- data[data20_indice, ]
head(data20)
m_20 <- stan_glm(y ~ x1 + x2_Treatment1 + x2_Treatment2 + 
                   x1:x2_Treatment1 + x1:x2_Treatment2, data = data20)
coef(m_20)
# (Intercept)    x1    x2_Treatment1  x2_Treatment2 x1:x2_Treatment1 x1:x2_Treatment2 
#-9.795209   6.206432    31.545891        -6.524821      1.479820        -0.917273 




# Parameter estimation across sample sizes
n_samples <- 10:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = length(coef(m1)))
colnames(bayes_coefficients) <- names(coef(m1))

for (i in n_samples) {
  sample_data <- data[1:i, ]
  m2_sample <- stan_glm(y ~ x1 + x2_Treatment1 + x2_Treatment2 + 
                          x1:x2_Treatment1 + x1:x2_Treatment2, data = sample_data)
  bayes_coefficients[i, ] <- coef(m2_sample)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_dummy_variable.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(ceiling(sqrt(ncol(bayes_coefficients))), ceiling(sqrt(ncol(bayes_coefficients))))) # Adjust layout for all parameters

for (j in seq_len(ncol(bayes_coefficients))) { # Use ncol() to iterate over all parameters
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2) # True parameter value
}
dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))





# Parameter estimation across sample sizes 10 rep
n_samples <- 20:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = length(coef(m1)))
colnames(bayes_coefficients) <- names(coef(m1))

for (i in n_samples) {
  
  sample_data <- data[1:i, ]
  draws <- replicate(10, {
    sample_indices <- sample(1:n, size = i, replace = FALSE) # Randomly select 'i' rows
    sample_data <- data[sample_indices, ]
    m_sample <- stan_glm(y ~ x1 + x2_Treatment1 + x2_Treatment2 + 
                           x1:x2_Treatment1 + x1:x2_Treatment2, data = sample_data)
    coef(m_sample)
  })
  bayes_coefficients[i, ] <- rowMeans(draws)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_dummy_variable_10rep.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(ceiling(sqrt(ncol(bayes_coefficients))), ceiling(sqrt(ncol(bayes_coefficients))))) # Adjust layout for all parameters

for (j in seq_len(ncol(bayes_coefficients))) { # Use ncol() to iterate over all parameters
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2) # True parameter value
}
dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))
