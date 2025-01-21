```{r, warning=FALSE, message=FALSE}
rm(list = ls()) #cleans your R space
options(stringsAsFactor = FALSE) #to stop R from automatically determining that your data is a factor

library(rstanarm)
library(ggplot2)
options(mc.cores = parallel::detectCores())

set.seed = 20690
```

# Homework 1
Isidora Silva-Valderrama
January 20th. 2025

## Part 1
**Simulate data for a linear regression problem (continuous predictor and response variables) of your choice and show a linear model can return the parameters you set. Your n should be 100 and you should set your error (sigma) to a level you would expect in an ecological or evolutionary biology study.**

Simulating predictive data for the linear regression: 
necrosis = alpha + beta1*diameter + error

y= alpha + beta1*x1 + error

```{r}
n <- 100 
x1 <- runif(n, min=0.5, max=3)#diamterer
alpha <- 0.5 #canes with smaller diameter than 0.5 are not considered 
beta1 <- 5
sigma <- 2
error <- rnorm(n,0,sigma)

y_pred = alpha + beta1*x1 + error
```
To see if we get the same estimates tha with our real data, we simulate some observed data: 
```{r}
y_obs <- y_pred + rnorm(n,0,sigma) #the values that already follow a linear pattern plus some noise
```
And we can "compare" them (visualize them): 
```{r}
plot(x1,y_pred, main = "Observed vs. Predicted", pch = 16, col = "blue", cex = 0.8) 
points(x1, y_obs,  pch = 17, col = "red", cex = 0.8)
legend("topleft", legend = c("Predicted (y_pred)", "Observed (y)"),
       col = c("blue", "red"), pch = c(16, 17), bty = "n")
```

```{r}
# Check if we can get the parameters back from the model with observed data
mod1 <- stan_glm(y_obs~x1)
summary(mod1)
```
The values of the parameters proposed / observed:
- alpha (intercept) 0.5 / **0.9**
- beta (slope) 5 / **5.1**
- sigma (stand. dev.) 2 / **3.0**

**Next, add a binary covariate (for example, a treatment and control) to your dataset that interacts with your other predictor variable. Make the effect size of this interaction HALF the effect size of the first predictor variable, and again test how well you can return the parameter values**
```{r}
#linear regression y= alpha + beta1*x1 + beta2*x2 + beta3*x1*x2 + error
#new covariate:
x2 <- rbinom(n, size = 1, prob = c(0.5,0.5)) #half with water, half without
beta2 <- 2
beta3 <- 2.5 #half beta1

#model
y_pred2 = alpha + beta1*x1 + beta2*x2 + beta3*x1*x2 + error 
#observed data generation:
y_obs2 <- y_pred2 + rnorm(n,0,sigma) 

#plot
par(mfrow = c(2, 2))
plot(y_pred2~x1*x2, main = "Predicted", pch = 16, col = "blue", cex = 0.8) 
plot(y_obs2~x1*x2, main = "Observed",  pch = 17, col = "red", cex = 0.8) 
```
```{r}
mod2 <- stan_glm(y_obs2~x1*x2)
summary(mod2)
```

The values of the parameters proposed / observed:
- alpha (intercept) 0.5 / **0.3**
- beta1 (slope) 5 / **5.2**
- beta2 (slope) 2 / **1.9**
- beta3 (slope inter.) 2.5 / **2.8**
- sigma (stand. dev.) 2 / **3.0**

## Part 2
**Step 2 -- Building on the above example, let's look at how the estimates of parameters returned from your linear model change with sample size, and make our analyses a little more robust. Using just 20% of your data, see how well you can return the parameters you set**

If I remove 80% of my observed data (n=100)
```{r}
# Randomly sample indices for the remaining 20% (which I am picking)
y_obs2_keep <- sample(seq_along(y_obs2), size = 20, replace = FALSE)

# Subset the vector to keep only the selected indices
y_obs2_2 <- y_obs2[y_obs2_keep]

#now, the model and estimates: 
x1 <- runif(20, min=0.5, max=3)
x2 <- rbinom(20, size = 1, prob = c(0.5,0.5))

mod3<-stan_glm(y_obs2_2 ~ x1*x2)
summary(mod3)
```

with only 20% of the data the differences between the proposed / observed data:
- alpha (intercept) 0.5 / **12.8**
- beta1 (slope) 5 / **0.6**
- beta2 (slope) 2 / **5.2**
- beta3 (slope inter.) 2.5 / **-2.3**
- sigma (stand. dev.) 2 / **7.7**

VERY different.

**Next, let's make a plot of sampling from your data from 1 to the n of your data (100) showing the estimated parameters for each sample size. Make a plot with your 1:100 on the horizontal axis and your estimated parameters on the vertical (you need either as many plots as parameters or a way to show all the different parameters on one plot). Compare how well the model does across the different parameters. Which is better or worse at and why?**

```{r}
#parameters
n <- c(5,10,20,40,30,50,60,70,80,100)
alpha <- 0.5
beta1 <- 5
beta2 <- 2
beta3 <- 2.5
sigma <- 2


#where to store the result:
results <- list()

# Loop over different sample sizes
for (i in n) {
  
  #the same simulations than above
  x1 <- runif(i, min=0.5, max=3)
  x2 <- rbinom(i, size = 1, prob = c(0.5,0.5))
  error <- rnorm(i, 0, sigma)
  y_pred4 = alpha + beta1*x1 + beta2*x2 + beta3*x1*x2 + error 
  y_obs4 <- y_pred4 + rnorm(i,0,sigma)
  
   # Check for variation in the data
  if (is.na(sd(x1)) || is.na(sd(y_obs4)) || is.na(sd(x2)) || sd(x1) == 0 || sd(y_obs4) == 0 || sd(x2) == 0) {
    next  # Skip this iteration if x1, x2 or y_obs4 has no variation or NA
  }
  
# Create a data rame with the first i observations
  data_subset <- data.frame(x1 = x1[1:i], x2 = x2[1:i], y_obs4 = y_obs4[1:i])
  # Fit the stan_glm() model
  fit <- stan_glm(y_obs4 ~ x1 * x2, data = data_subset)
  # Extract summary statistics
  fit_summary <- as.data.frame(summary(fit))
  # Filter estimates for specific parameters
  estimates <- fit_summary[rownames(fit_summary) %in% c("(Intercept)", "x1", "x2", "x1:x2", "sigma"), 
                           "mean"]
  # Add a column for sample size
  estimates$sample_size <- i
  # Add parameter names as a column
  estimates$parameter <- rownames(estimates)
  # Append the results to the list
  results[[i]] <- estimates
}

# Combine all results into a single data frame
final_results <- do.call(rbind, results)
colnames(final_results) <- c("Intercept", "beta1", "beta2", "beta3", "sigma", "sample_size")
final_results <- as.data.frame(final_results)
final_results[] <- lapply(final_results, as.numeric)


#Plot
ggplot(final_results, aes(x=sample_size, y=Intercept)) + 
  geom_line()+
  annotate("label", x = 105, y = -1.0, label = "Intercept", size= 3) +
  geom_line(aes(y=beta1), color="orange")+
  annotate("label", x = 105, y = 6, label = "beta1", size= 3, color="orange") +
  geom_line(aes(y=beta2), color="red")+
  annotate("label", x = 105, y = 4.3, label = "beta2", size= 3, color="red") +
  geom_line(aes(y=beta3), color="pink")+
  annotate("label", x = 105, y = 1.0, label = "beta3", size= 3, color="pink") +
  geom_line(aes(y=sigma), color="purple")+
  annotate("label", x = 105, y = 2.64, label = "sigma", size= 3, color="purple") +
  theme_classic()+
  ylim(-10,13)+
  labs(y= "Parameter", x="sample size")
```

Around sample 35-40, the parameters show less variation with sample size. The intercept (alpha), and beta1 show the highest variation with small sample size, and it is produced by the variation or noise included with each sample (there could be outliers too). Therefore, the estimates are more "sensitive" to the random noise and the variation from the sampling process or error. Sigma is the parameter that stays the most consistent.  

**Let's improve our sampling now. So far we have take just ONE draw from our set of parameters which means our sigma term has some Monte Carlo error in it, so let's take 10 draws each time (so we need to set up a loop or such in R that samples from 1:n and at each step it repeats simulating the data and getting the estimates 10 times). Re-make your plot. (If you get stuck here for a while, don't panic, but move onto Steps 3-4.)**
```{r}
#using the same parameters
#where to store the result:
results <- list()

# Loop over different sample sizes
for (i in n) {
  
   # temporary list to store results for 10 draws at each step
  temp_results <- list()
  
  # Nested mini loop to get the 10 draws each time
  for (j in 1:10) {
    # Simulate data
    x1 <- runif(i, min = 0.5, max = 3)
    x2 <- rbinom(i, size = 1, prob = 0.5)
    error <- rnorm(i, 0, sigma)
    y_pred4 <- alpha + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2 + error
    y_obs4 <- y_pred4 + rnorm(i, 0, sigma)
   # Check for variation in the data
  if (is.na(sd(x1)) || is.na(sd(y_obs4)) || is.na(sd(x2)) || sd(x1) == 0 || sd(y_obs4) == 0 || sd(x2) == 0) {
    next  # Skip this iteration if x1, x2 or y_obs4 has no variation or NA
  }
  
# Create a data rame with the first i observations
  data_subset <- data.frame(x1 = x1[1:i], x2 = x2[1:i], y_obs4 = y_obs4[1:i])
  # Fit the stan_glm() model
  fit <- stan_glm(y_obs4 ~ x1 * x2, data = data_subset)
  # Extract summary statistics
  fit_summary <- as.data.frame(summary(fit))
  # Filter estimates for specific parameters
  estimates <- fit_summary[rownames(fit_summary) %in% c("(Intercept)", "x1", "x2", "x1:x2", "sigma"), 
                           "mean"]
  # Add a column for sample size
  estimates$sample_size <- i
  # Add parameter names as a column
  estimates$parameter <- rownames(estimates)
  # Append the results to the temp. list
  temp_results[[j]] <- estimates
}

    # Combine results for the 10 draws at this sample size
  if (length(temp_results) > 0) {
    results[[i]] <- do.call(rbind, temp_results)
  }
}

# Combine all results into the final dataframe
final_results <- do.call(rbind, results)
colnames(final_results) <- c("Intercept", "beta1", "beta2", "beta3", "sigma", "sample_size")
final_results <- as.data.frame(final_results)
final_results[] <- lapply(final_results, as.numeric)


#Plot
ggplot(final_results, aes(x=sample_size, y=Intercept)) + 
  geom_jitter()+
  annotate("label", x = 110, y = -1.4, label = "Intercept", size= 3) +
  geom_jitter(aes(y=beta1), color="orange")+
  annotate("label", x = 110, y = 6, label = "beta1", size= 3, color="orange") +
  geom_jitter(aes(y=beta2), color="red")+
  annotate("label", x = 110, y = 4.3, label = "beta2", size= 3, color="red") +
  geom_jitter(aes(y=beta3), color="pink")+
  annotate("label", x = 110, y = 1.0, label = "beta3", size= 3, color="pink") +
  geom_jitter(aes(y=sigma), color="purple")+
  annotate("label", x = 110, y = 2.64, label = "sigma", size= 3, color="purple") +
  theme_classic()+
  labs(y= "Parameter", x="sample_size")
```

It definitely shows less variation after sample size ~40.

## Part 3
**Now, repeat the above but change your error term by 50%. Compare these results to what you found before and explain why you think they changed**
```{r}
#parameters
n <- c(5,10,20,40,30,50,60,70,80,100)
alpha <- 0.5
beta1 <- 5
beta2 <- 2
beta3 <- 2.5
sigma <- 4 #sigma, but increased 50% (was 2)

#where to store the result:
results <- list()

# Loop over different sample sizes
for (i in n) {
  
   # temporary list to store results for 10 draws at each step
  temp_results <- list()
  
  # Nested mini loop to get the 10 draws each time
  for (j in 1:10) {
    # Simulate data
    x1 <- runif(i, min = 0.5, max = 3)
    x2 <- rbinom(i, size = 1, prob = 0.5)
    error <- rnorm(i, 0, sigma)
    y_pred4 <- alpha + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2 + error
    y_obs4 <- y_pred4 + rnorm(i, 0, sigma)
   # Check for variation in the data
  if (is.na(sd(x1)) || is.na(sd(y_obs4)) || is.na(sd(x2)) || sd(x1) == 0 || sd(y_obs4) == 0 || sd(x2) == 0) {
    next  # Skip this iteration if x1, x2 or y_obs4 has no variation or NA
  }
  
# Create a data rame with the first i observations
  data_subset <- data.frame(x1 = x1[1:i], x2 = x2[1:i], y_obs4 = y_obs4[1:i])
  # Fit the stan_glm() model
  fit <- stan_glm(y_obs4 ~ x1 * x2, data = data_subset)
  # Extract summary statistics
  fit_summary <- as.data.frame(summary(fit))
  # Filter estimates for specific parameters
  estimates <- fit_summary[rownames(fit_summary) %in% c("(Intercept)", "x1", "x2", "x1:x2", "sigma"), 
                           "mean"]
  # Add a column for sample size
  estimates$sample_size <- i
  # Add parameter names as a column
  estimates$parameter <- rownames(estimates)
  # Append the results to the temp. list
  temp_results[[j]] <- estimates
}

    # Combine results for the 10 draws at this sample size
  if (length(temp_results) > 0) {
    results[[i]] <- do.call(rbind, temp_results)
  }
}

# Combine all results into the final dataframe
final_results <- do.call(rbind, results)
colnames(final_results) <- c("Intercept", "beta1", "beta2", "beta3", "sigma", "sample_size")
final_results <- as.data.frame(final_results)
final_results[] <- lapply(final_results, as.numeric)


#Plot
ggplot(final_results, aes(x=sample_size, y=Intercept)) + 
  geom_jitter()+
  annotate("label", x = 110, y = -1.4, label = "Intercept", size= 3) +
  geom_jitter(aes(y=beta1), color="orange")+
  annotate("label", x = 110, y = 6, label = "beta1", size= 3, color="orange") +
  geom_jitter(aes(y=beta2), color="red")+
  annotate("label", x = 110, y = 4.3, label = "beta2", size= 3, color="red") +
  geom_jitter(aes(y=beta3), color="pink")+
  annotate("label", x = 110, y = 1.0, label = "beta3", size= 3, color="pink") +
  geom_jitter(aes(y=sigma), color="purple")+
  annotate("label", x = 110, y = 2.64, label = "sigma", size= 3, color="purple") +
  theme_classic()+
  labs(y= "Parameter", x="sample_size")

```

The variation of the parameters (specially beta2) with a higher sigma, increased drastically in the models with low sample size. Therefore, if the error or the noise in the data is high, a higher sample size is needed to stabilize the parameters.

## Part 4
**Review the sample datasets and pick one you want to use to go through the complete workflow with (which is the homework for next week). Ideally, you'll get to work on that one, but I may ask people to pick a new one so everyone is not using the same one.**

I'll like to work with the carnivore dataset as first preference. If not, then the damselflies. 
