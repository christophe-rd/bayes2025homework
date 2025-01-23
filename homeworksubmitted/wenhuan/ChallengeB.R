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

# Fit a Bayesian model including all interactions
m1 <- stan_glm(y ~ x1 + treatment1_Treatment + treatment2_Treatment + 
                 x1:treatment1_Treatment + x1:treatment2_Treatment + 
                 treatment1_Treatment:treatment2_Treatment + 
                 x1:treatment1_Treatment:treatment2_Treatment, 
               data = data)
summary(m1)
coef(m1)

# Repeat for a smaller sample size (20%)
set.seed(123)
data20_indices <- sample(1:n, size = 20, replace = FALSE)
data20 <- data[data20_indices, ]
head(data20)
m_20 <- stan_glm(y ~ x1 + treatment1_Treatment + treatment2_Treatment + 
                   x1:treatment1_Treatment + x1:treatment2_Treatment + 
                   treatment1_Treatment:treatment2_Treatment + 
                   x1:treatment1_Treatment:treatment2_Treatment, 
                 data = data20)
coef(m_20)

# Parameter estimation across sample sizes
n_samples <- 20:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = length(coef(m1)))
colnames(bayes_coefficients) <- names(coef(m1))

for (i in n_samples) {
  sample_data <- data[1:i, ]
  m2_sample <- stan_glm(y ~ x1 + treatment1_Treatment + treatment2_Treatment + 
                          x1:treatment1_Treatment + x1:treatment2_Treatment + 
                          treatment1_Treatment:treatment2_Treatment + 
                          x1:treatment1_Treatment:treatment2_Treatment, 
                        data = sample_data)
  bayes_coefficients[i, ] <- coef(m2_sample)
}

# Plot estimated parameters across sample sizes
tiff(filename = "parameter_estimates_full_factorial.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(ceiling(sqrt(ncol(bayes_coefficients))), ceiling(sqrt(ncol(bayes_coefficients)))))

for (j in seq_len(ncol(bayes_coefficients))) {
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2)
}
dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))

# Parameter estimation with 10 repetitions per sample size
n_samples <- 20:n
bayes_coefficients <- matrix(NA, nrow = n, ncol = length(coef(m1)))
colnames(bayes_coefficients) <- names(coef(m1))

for (i in n_samples) {
  draws <- replicate(10, {
    sample_indices <- sample(1:n, size = i, replace = FALSE)
    sample_data <- data[sample_indices, ]
    m_sample <- stan_glm(y ~ x1 + treatment1_Treatment + treatment2_Treatment + 
                           x1:treatment1_Treatment + x1:treatment2_Treatment + 
                           treatment1_Treatment:treatment2_Treatment + 
                           x1:treatment1_Treatment:treatment2_Treatment, 
                         data = sample_data)
    coef(m_sample)
  })
  bayes_coefficients[i, ] <- rowMeans(draws)
}

# Plot estimated parameters across sample sizes with 10 repetitions
tiff(filename = "parameter_estimates_full_factorial_10rep.tiff", width = 10, height = 7, units = "in", res = 300)
par(mfrow = c(ceiling(sqrt(ncol(bayes_coefficients))), ceiling(sqrt(ncol(bayes_coefficients)))))

for (j in seq_len(ncol(bayes_coefficients))) {
  plot(n_samples, na.omit(bayes_coefficients[, j]), type = "l", 
       main = paste("Parameter:", colnames(bayes_coefficients)[j]),
       xlab = "Sample Size", ylab = "Estimated Value",
       col = "blue")
  abline(h = coef(m1)[j], col = "red", lwd = 2, lty = 2)
}
dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))


