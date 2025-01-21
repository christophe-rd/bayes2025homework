### FRST507 - Homework 1
## Hannah Bates - Jan 2025

# load stan pkg
library(rstanarm)

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)


### STEP 1 ----
## SIMULATING DATA - no treatments
n <- 100 #sample size
a <- 1 # intercept
b <- 2 # slope; "effect size"
sigma <- 40 #variance

# y = a + bx + error
x <- runif(n,0,100) # draw random x values from uniform dist btw 0 and 100
ypred <- a + b*x # predicted values
yobs <- ypred+rnorm(n,0,sigma) # observed values - including error
plot(x, yobs)

# Fit linear model using stan
m1 <- stan_glm(yobs~x)
summary(m1)
abline(a=coef(m1)[1], b=coef(m1)[2], col = "red")

# Results:
# The intercept from the model (8.1) is not very close to the given parameter
# (1), however the slope given by the model is exact (2.0)  and the variance 
# is very similar (41.2).


## SIMULATING DATA - w binary covariate: treatment (1) and control (0)
rm(list=ls())

n <- 100 #sample size
a <- 1 # intercept
b1 <- 2 # slope; "effect size"
b2 <- -2 # slope for treatment effect
b3 <- 1 # effect size of interaction
sigma <- 40 #variance

x1 <- runif(n,0,100) # draw random x values from uniform dist btw 0 and 100
x2 <- rbinom(n,1,0.5) # 0 for control and 1 for treatment: split 50/50 between tanks
ypred <- a + b1*x1 + b2*x2 + b3*x1*x2 # predicted values
yobs <- ypred+rnorm(n,0,sigma) # observed values - including error
plot(x1, yobs)

# Fit linear model using stan
m2 <- stan_glm(yobs ~ x1*x2) # only need to write interaction; will include indv effects
summary(m2)

# Results:
# The intercept from the model (10.6) is not close to the given 
# parameter (1). The slope of the first predictor variable x1 (2), as well as the
# interaction effect size (1), are very close to the values given by the model
# (1.9 and 1.2, respectively). The slope for the treatment effect (-2) was not 
# accurately returned by the model, which estimated -21.3. Sigma was close at 44.2.


### STEP 2 ----
## CHANGING SAMPLE SIZE

# Create subset of 20%:
indicies <- sample(1:n, 0.2*n, replace = FALSE) # sampling 20% of original sample size
x1subset <- x1[indicies] # creating subsets from randomly sampled indicies
x2subset <- x2[indicies]
yobssubset <- yobs[indicies]
plot(x1subset, yobssubset)

# Fit linear model using stan
m3 <- stan_glm(yobssubset ~ x1subset*x2subset) # only need to write interaction; will include indv effects
summary(m3)

# Results:
# As previously, the model returned similar estimates for the x1 slope, 
# the interaction effect, and sigma: however, these values were not as accurate 
# as with the full sample size. The values for intercept and treatment slope which
# were returned were very different from the initial parameters.

## PLOT CHANGING SAMPLE SIZE

# Create empty dataframe to store coefs to plot:
coefs <- data.frame(matrix(ncol = 5, nrow = n))
colnames(coefs) <- c("intercept", "x1slope", "x2slope", "interaction", "sigma")

# Create for loop:
for (i in seq(10, n, by = 10)) { # increasing in intervals of 10 to limit computing time
  indicies2 <- sample(1:n, i, replace = FALSE) # sampling i values from original sample size
  x1subset <- x1[indicies2] # creating subsets from randomly sampled indicies
  x2subset <- x2[indicies2]
  yobssubset <- yobs[indicies2]
  m4 <- stan_glm(yobssubset ~ x1subset*x2subset) # fitting model at each sample size
  coefs$intercept[i] <- coef(m4)[1] # storing coefficient values
  coefs$x1slope[i] <- coef(m4)[2]
  coefs$x2slope[i] <- coef(m4)[3]
  coefs$interaction[i] <- coef(m4)[4]
  coefs$sigma[i] <- summary(m4)["sigma","mean"]
}

# Plot results
par(mfrow = c(2, 3))
plot(1:n, coefs$intercept)
plot(1:n, coefs$x1slope)
plot(1:n, coefs$x2slope)
plot(1:n, coefs$interaction)
plot(1:n, coefs$sigma)

# Results/Comparison:
# The model appears to estimate the first covariate and the interaction parameters 
# best at lower sample sizes, compared with the parameters for the intercept and
# the treatment covariate. The latter two have very wide ranges at smaller sample sizes,
# whereas the former two begin with much smaller ranges (between 0-3). It appears that
# the interaction term becomes more accurate more quickly as sample size increases.
# The estimated sigma is not very accurate although settles after a sample size greater
# than the real sigma.


## MULTIPLE SAMPLING

# Empty data frames for each parameter
intercept_df <- data.frame(matrix(ncol = 10, nrow = n))
x1slope_df <- data.frame(matrix(ncol = 10, nrow = n))
x2slope_df <- data.frame(matrix(ncol = 10, nrow = n))
interaction_df <- data.frame(matrix(ncol = 10, nrow = n))
sigma_df <- data.frame(matrix(ncol = 10, nrow = n))

# Create nested for loop
for (i in seq(10, n, by = 10)) { # increasing in intervals of 10 to limit computing time
  for (j in 1:10) { # simulating data 10 times at each sample size
    x1 <- runif(i,0,100)
    x2 <- rbinom(i,1,0.5)
    ypred <- a + b1*x1 + b2*x2 + b3*x1*x2
    yobs <- ypred+rnorm(i,0,sigma) # new yobs with each iteration
    m5 <- stan_glm(yobs ~ x1*x2) 
    intercept_df[i,j] <- coef(m5)[1] # storing coefficient values for each size and simulation
    x1slope_df[i,j] <- coef(m5)[2]
    x2slope_df[i,j] <- coef(m5)[3]
    interaction_df[i,j] <- coef(m5)[4]
    sigma_df[i,j] <- summary(m5)["sigma","mean"]
  }
}

# Plot results
par(mfrow = c(2, 3))
matplot(intercept_df, type = "p", pch = 1) # different colours for different simulated data
matplot(x1slope_df, type = "p", pch = 1)
matplot(x2slope_df, type = "p", pch = 1)
matplot(interaction_df, type = "p", pch = 1)
matplot(sigma_df, type = "p", pch = 1)


### STEP 3 ----
## CHANGE ERROR BY 50%

# Increase sigma by 50%
sigma <- 60 #variance - originally: 40

# Repeat simulating data with new sigma value:
x1 <- runif(n,0,100) # step 2 changed my x1/x2 values
x2 <- rbinom(n,1,0.5)
ypred <- a + b1*x1 + b2*x2 + b3*x1*x2
yobs <- ypred+rnorm(i,0,sigma) # new sigma!

## Repeat: PLOT CHANGING SAMPLE SIZE

# Create empty dataframe to store coefs to plot:
coefs2 <- data.frame(matrix(ncol = 5, nrow = n))
colnames(coefs2) <- c("intercept", "x1slope", "x2slope", "interaction", "sigma")

# Create for loop:
for (i in seq(10, n, by = 10)) {
  indicies3 <- sample(1:n, i, replace = FALSE) # sampling i values
  x1subset <- x1[indicies3] # creating subsets from randomly sampled indicies
  x2subset <- x2[indicies3]
  yobssubset <- yobs[indicies3]
  m6 <- stan_glm(yobssubset ~ x1subset*x2subset)
  coefs2$intercept[i] <- coef(m6)[1]
  coefs2$x1slope[i] <- coef(m6)[2]
  coefs2$x2slope[i] <- coef(m6)[3]
  coefs2$interaction[i] <- coef(m6)[4]
  coefs2$sigma[i] <- summary(m6)["sigma","mean"]
}

# Plot results
par(mfrow = c(2, 3))
plot(1:n, coefs2$intercept)
plot(1:n, coefs2$x1slope)
plot(1:n, coefs2$x2slope)
plot(1:n, coefs2$interaction)
plot(1:n, coefs2$sigma)

# Results: 
# With increased error, the model did not generate as similar values as in Step 2.
# As previously, the model managed to estimate the first covariate and the 
# interaction parameters more accurately than other parameters at smaller sample
# sizes. The intercept was most difficult to estimate at smaller sample sizes,
# and the sigma value managed to be reasonably estimated once the sample size
# grew larger than the true sigma value. I think these results changed as more error
# was put into the system, because the true trends became more masked and less 
# easy to tease out at smaller sample sizes.


### CHALLENGE C - Extract & Plot Uncertainty ----

# Empty data frames for each parameter
intercept_df <- data.frame(matrix(ncol = 3, nrow = n))
x1slope_df <- data.frame(matrix(ncol = 3, nrow = n))
x2slope_df <- data.frame(matrix(ncol = 3, nrow = n))
interaction_df <- data.frame(matrix(ncol = 3, nrow = n))
sigma_df <- data.frame(matrix(ncol = 3, nrow = n))

# Create for loop:
i = 10 # not sure if this is actually necessary but just in case
for (i in seq(10, n, by = 10)) {
  indicies3 <- sample(1:n, i, replace = FALSE) # sampling i values from original sample size
  x1subset <- x1[indicies3] # creating subsets from randomly sampled indicies
  x2subset <- x2[indicies3]
  yobssubset <- yobs[indicies3]
  m7 <- stan_glm(yobssubset ~ x1subset*x2subset)
  intercept_df[i,1] <- coef(m7)[1] # storing the means
  x1slope_df[i,1] <- coef(m7)[2]
  x2slope_df[i,1] <- coef(m7)[3]
  interaction_df[i,1] <- coef(m7)[4]
  sigma_df[i,1] <- summary(m7)["sigma","mean"]
  intercept_df[i,2] <- posterior_interval(m7, prob = 0.5)[1,1] # storing the 25% uncertainty
  x1slope_df[i,2] <- posterior_interval(m7, prob = 0.5)[2,1]
  x2slope_df[i,2] <- posterior_interval(m7, prob = 0.5)[3,1]
  interaction_df[i,2] <- posterior_interval(m7, prob = 0.5)[4,1]
  sigma_df[i,2] <- posterior_interval(m7, prob = 0.5)[5,1]
  intercept_df[i,3] <- posterior_interval(m7, prob = 0.5)[1,2] # storing the 75% uncertainty
  x1slope_df[i,3] <- posterior_interval(m7, prob = 0.5)[2,2]
  x2slope_df[i,3] <- posterior_interval(m7, prob = 0.5)[3,2]
  interaction_df[i,3] <- posterior_interval(m7, prob = 0.5)[4,2]
  sigma_df[i,3] <- posterior_interval(m7, prob = 0.5)[5,2]
}

# Plot means with error bars
par(mfrow = c(2, 3))
plot(1:n, intercept_df[,1])
arrows(1:n, intercept_df[,2], 1:n, intercept_df[,3], angle = 90, code = 3, length = 0.1)
plot(1:n, x1slope_df[,1])
arrows(1:n, x1slope_df[,2], 1:n, x1slope_df[,3], angle = 90, code = 3, length = 0.1)
plot(1:n, x2slope_df[,1])
arrows(1:n, x2slope_df[,2], 1:n, x2slope_df[,3], angle = 90, code = 3, length = 0.1)
plot(1:n, interaction_df[,1])
arrows(1:n, interaction_df[,2], 1:n, interaction_df[,3], angle = 90, code = 3, length = 0.1)
plot(1:n, sigma_df[,1])
arrows(1:n, sigma_df[,2], 1:n, sigma_df[,3], angle = 90, code = 3, length = 0.1)

# Comparing across parameters at n = 100
plot(1:5, summary(m7)[1:5,"mean"], xlab = "Parameters: a, b1, b2, b3, sig", ylab = "mean parameter estimates", ylim = c(-20,70))
arrows(1:5, posterior_interval(m7, prob = 0.5)[,1], 1:5, posterior_interval(m7, prob = 0.5)[,2], angle = 90, code = 3, length = 0.1)

# The above plot shows how little error the first covariate and the interaction 
# parameter have in comparison to the other parameters, even at the maximum sample size.

### STEP 4 ----

# I'd love to work with the carnivore data!