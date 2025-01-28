## Started 27 January 2025
## By Ken

## Bayes 2025, homework 2, workflow attempt, only steps 1, 2, and 4

rm(list = ls())
options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
options(mc.cores = parallel::detectCores())

# set working directory
if(length(grep("Ken", getwd()) > 0)){
  setwd("/Users/Ken Michiko Samson/Documents/2024W2/Bayes2025/bayes2025homework/homeworksubmitted2/ken")
}

## Step 1

# format data
d <- read.delim("input/fastswim.dta", sep = "")
d$sp <- NA
for(i in 1:length(unique(d$species))){
  d$sp[which(d$species == unique(d$species)[i])] <- i
}

# data visualization
par(mfrow = c(3, 3))
plot(d$wetmass, d$speed, col = d$sp)
plot(d$swimbeat, d$speed)
plot(d$tailamp, d$speed)

plot(d$ldh, d$speed)
plot(d$pk, d$speed)
plot(d$ak, d$speed)

plot(d$area, d$speed)
plot(d$perim, d$speed)
plot(d$majaxis, d$speed)

plot(d$breadth, d$speed, col = d$sp)

par(mfrow = c(1, 1))

d$wetmass <- as.numeric(d$wetmass)
d$speed <- as.numeric(d$speed)
d$sp <- as.numeric(d$sp)

# model speed as a logistic function of wetmass
# max / (1 + e^(-k(x - midpoint)))

## Step 2

# logistic model paramters
N_sp <- 11

# Population parameters
mmax <- 30
smax <- 3 / 2.32
mk <- 0.2
sk <- 0.05 / 2.32
mmp <- 13
smp <- 4 / 2.32
error <- 5

# Define species-level parameters
maxspd <- rnorm(N_sp, mmax, smax)
k <- rnorm(N_sp, mk, sk)
mp <- rnorm(N_sp, mmp, smp)

# Simulate data
N <- 10000

wm <- rep(NA, N)
sp <- round(runif(N, 0.5, N_sp + 0.5))

lp <- rnorm(N_sp, 2.5, 0.5 / 2.32)
sd <- rnorm(N_sp, 1, 0.5 / 2.32)

for(i in 1:N){
  wm[i] <- rlnorm(1, lp[sp[i]], sd[sp[i]])
}

spd <- maxspd[sp] / (1 + exp(-k[sp] * (wm - mp[sp])))
spd_obs <- abs(rnorm(N, spd, error))
plot(wm, spd_obs, col = sp, xlim = c(0, 50), ylim = c(0, 40))

data <- list("N" = N, "N_sp" = N_sp, "wm" = wm, "sp" = sp, "spd" = spd_obs)

fit <- stan(file = "stan/logistic_hier.stan", data = data,
            warmup = 2000, iter = 3000, refresh = 100, chains = 4)
saveRDS(fit, "simu_data_fit.rds")

#launch_shinystan(fit)
# checking the model with simulated fit yielded Rhat = 1, n_eff > 3500 (87.5%), and no divergences
# population parameters were returned with good approximations
# however, divergences were observed with N = 100 and 1000

## Step 4

dwv <- subset(d, !is.na(sp) & !is.na(wetmass) & !is.na(speed))

data <- list("N" = nrow(dwv), "N_sp" = length(unique(dwv$sp)), "wm" = dwv$wetmass, "sp" = dwv$sp, "spd" = dwv$speed)

fit <- stan(file = "stan/logistic_hier.stan", data = data,
            warmup = 2000, iter = 3000, refresh = 100, chains = 4)
saveRDS(fit, "fit.rds")

#launch_shinystan(fit)
# checking the model with simulated fit yielded Rhat = 1, n_eff > 1000 (25%), and no divergences
# significant sample correlation were not observed
