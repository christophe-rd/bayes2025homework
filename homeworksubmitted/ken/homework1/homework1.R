## Started 19 January 2025
## By Ken

## Bayes 2025, homework 1, simulating data

rm(list = ls())
options(stringsAsFactors = FALSE)

library(rstanarm)
options(mc.cores = parallel::detectCores())

# step 1

# parameters for linear regression with one predictor
a <- 50 #intercept
b <- 5 #effect size
e <- 50 #error

# sample size
N <- 100

# simulate data
x <- runif(N, 0, 100)
y <- rnorm(N, a + b*x, e)

m <- stan_glm(y ~ x)

plot(x, y)
abline(m$coefficients[1], m$coefficients[2])
print(paste0("Intercept (a): ", a, "; Estimate: ", m$coefficients[1]))
print(paste0("Slope (b): ", b, "; Estimate: ", m$coefficients[2]))
print(paste0("Error (e): ", e, "; Estimate: ", summary(m)[3, 1]))

# intercept fluctuates +- 10
# slope fluctuates +- 0.5
# error fluctuates +- 5


# parameters for linear regression with two predictors: one continuous (b1) and one binary (b2)
# model:
# y = a + b1*x1 + b2*x2 + b3*x1*x2
# e ~ N(0, e)
a <- 50 #intercept
b1 <- 5 #effect size of first variable
b2 <- 20 #effect size of second variable, i call it added intercept
b3 <- 2.5 #interaction
e <- 50 #error

# simulate data
x1 <- runif(N, 0, 100)
x2 <- rbinom(N, 1, 0.5)
y <- rnorm(N, a + b1*x1 + b2*x2 + b3*x1*x2, e)

m <- stan_glm(y ~ x1*x2)
m$coefficients

plot(x1, y)
abline(m$coefficients[1], m$coefficients[2], col = 'red')
abline(m$coefficients[1] + m$coefficients[3],
       m$coefficients[2] + m$coefficients[4], col = 'blue')
print(paste0("Intercept (a): ", a, "; Estimate: ", m$coefficients[1]))
print(paste0("Slope (b1): ", b1, "; Estimate: ", m$coefficients[2]))
print(paste0("Added intercept (b2): ", b2, "; Estimate: ", m$coefficients[3]))
print(paste0("Interaction (b3): ", b3, "; Estimate: ", m$coefficients[4]))
print(paste0("Error (e): ", e, "; Estimate: ", summary(m)[5, 1]))

# intercept fluctuates +- 20, more than before 
# slope fluctuates +- 0.5 still
# added intercept fluctuates wildly, usually with intercept estimates as well
# interaction fluctuates +- 1
# error fluctuates +- 5 still

# step 2

# loop to vary sample size
params <- array(NA, c(5, 100))
for(i in 2:100){
  print(i)
  N <- i
  
  x1 <- runif(N, 0, 100)
  x2 <- rbinom(N, 1, 0.5)
  y <- rnorm(N, a + b1*x1 + b2*x2 + b3*x1*x2, e)
  
  m <- stan_glm(y ~ x1*x2)
  
  params[1:4, i] <- m$coefficients
  params[5, i] <- summary(m)[5, 1]
}

par(mfrow = c(2, 3))
plot(1:100, params[1, ], ylim = c(10, 90),
     xlab = "Sample Size (N)", ylab = "Intercept (a)")
abline(a = a, b = 0)
plot(1:100, params[2, ], ylim = c(0, 10),
     xlab = "Sample Size (N)", ylab = "Effect Size (b1)")
abline(a = b1, b = 0)
plot(1:100, params[3, ], ylim = c(-20, 60),
     xlab = "Sample Size (N)", ylab = "Added intercept (b2)")
abline(a = b2, b = 0)
plot(1:100, params[4, ], ylim = c(0, 5),
     xlab = "Sample Size (N)", ylab = "Interaction (b3)")
abline(a = b3, b = 0)
plot(1:100, params[5, ], ylim = c(30, 70),
     xlab = "Sample Size (N)", ylab = "Error (e)")
abline(a = e, b = 0)

# i can't fit the model with only one data point
# it was left as NA

# only the effect size of the first variable seems to congregate to some range of the true value
# the interaction has less of a converging behaviour while the rest seem to be quite scattered

# this might be because the error allows for most of the data points in each group
# (binary covariate 0 and 1) to be around the same values for low values of x
# which may lead to difficulty in separating the interaction with the effect size of the
# continuous variable, essential separating the two groups

# loop to vary sample size with 10 estimates for each x
params <- array(NA, c(5, 100))
for(i in 2:100){
  N <- i
  x1 <- runif(N, 0, 100)
  x2 <- rbinom(N, 1, 0.5)
  
  x1 <- rep(x1, each = 10)
  x2 <- rep(x2, each = 10)
  
  y <- rnorm(N*10, a + b1*x1 + b2*x2 + b3*x1*x2, e)
  
  m <- stan_glm(y ~ x1*x2)
  
  params[1:4, i] <- m$coefficients
  params[5, i] <- summary(m)[5, 1]
}

par(mfrow = c(2, 3))
plot(1:100, params[1, ], ylim = c(10, 90),
     xlab = "Sample Size (N)", ylab = "Intercept (a)")
abline(a = a, b = 0)
plot(1:100, params[2, ], ylim = c(0, 10),
     xlab = "Sample Size (N)", ylab = "Effect Size (b1)")
abline(a = b1, b = 0)
plot(1:100, params[3, ], ylim = c(-20, 60),
     xlab = "Sample Size (N)", ylab = "Added intercept (b2)")
abline(a = b2, b = 0)
plot(1:100, params[4, ], ylim = c(0, 5),
     xlab = "Sample Size (N)", ylab = "Interaction (b3)")
abline(a = b3, b = 0)
plot(1:100, params[5, ], ylim = c(30, 70),
     xlab = "Sample Size (N)", ylab = "Error (e)")
abline(a = e, b = 0)

# estimates for all 5 parameters improved

# step 3
e <- 25

# loop to vary sample size with 10 estimates for each x
params <- array(NA, c(5, 100))
for(i in 2:100){
  N <- i
  x1 <- runif(N, 0, 100)
  x2 <- rbinom(N, 1, 0.5)
  
  x1 <- rep(x1, each = 10)
  x2 <- rep(x2, each = 10)
  
  y <- rnorm(N*10, a + b1*x1 + b2*x2 + b3*x1*x2, e)
  
  m <- stan_glm(y ~ x1*x2)
  
  params[1:4, i] <- m$coefficients
  params[5, i] <- summary(m)[5, 1]
}

par(mfrow = c(2, 3))
plot(1:100, params[1, ], ylim = c(10, 90),
     xlab = "Sample Size (N)", ylab = "Intercept (a)")
abline(a = a, b = 0)
plot(1:100, params[2, ], ylim = c(0, 10),
     xlab = "Sample Size (N)", ylab = "Effect Size (b1)")
abline(a = b1, b = 0)
plot(1:100, params[3, ], ylim = c(-20, 60),
     xlab = "Sample Size (N)", ylab = "Added intercept (b2)")
abline(a = b2, b = 0)
plot(1:100, params[4, ], ylim = c(0, 5),
     xlab = "Sample Size (N)", ylab = "Interaction (b3)")
abline(a = b3, b = 0)
plot(1:100, params[5, ], ylim = c(0, 50),
     xlab = "Sample Size (N)", ylab = "Error (e)")
abline(a = e, b = 0)

# the parameter estimates converge to a smaller range of the true value of the parameters
# this ties into how the error might make the distribution of values of the two groups overlap
# with each other

# by decreasing the error, the overlap is lessened and the effect size of the continuous variable
# and the interaction are estimated better, and by extension the intercepts (intercept and effect size
# of binary variable)

# step 4

# the dataset of damselfly swimming performance seems to be interesting
