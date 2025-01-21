#Bayesian Course -- Homework 1
#By Jenna Melanson
#Prof Lizzie Wolkovich

#README
# All functions for simulation and monte carlo are stored in src/simulationfns.R
# All figures are saved in figures folder

#housekeeping
rm(list=ls())

#set your working directory here!! should be the folder "jenna" in order for
#function sourcing and figure writing to work properly
setwd("/Users/jenna1/Documents/UBC/Coursework/bayesian_wolkovich")
options(stringsAsFactors = FALSE)
source("src/simulationfns.R")

#load packages
library(ggplot2)
library(rstanarm)
library(dplyr)
library(gridExtra)


#####################################
#### PART 1
####################################

#Simulate data for a linear regression problem and show a linear model can return the 
#parameters you set. n = 100, biologically realistic sigma

# Use the function "simulateLinRegData"
n = 100
simOnePredictor = simulateLinRegData(sample_size = n,
                                     intercept = 10,
                                     coefs = c(2),
                                     predictors = list(runif(n, 0, 1)),
                                     sigma = 0.7)

#plot to see if it looks *real*
ggplot(simOnePredictor, aes(x = pred1, y = yobs)) +
  geom_point()

#sigma = 0.7 gives an okay spread -- still pretty *nice* for ecological data imo

#fit model using rstanarm
fit = rstanarm::stan_glm(yobs ~ pred1, data = simOnePredictor)

#plot with fit line
ggplot(simOnePredictor, aes(x = pred1, y = yobs)) +
  geom_point() +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], col = 'red')
#looks pretty good

#Model coefficients are stored here:
coef(fit)

#Based on my simulation, the estimate intercept is 9.9 (true intercept = 10)
#and the estimated slope is 2.1 (true = 2)


#Next, add a binary covariate to your dataset that interacts with your other 
#predictor variable, effect size of interaction = HALF the effect size of pred1
b1 = 2
b2 = 0
interaction_effect = b1/2

#simulate data
n = 100
simTwoPredInteraction = simulateLinRegData(sample_size = n,
                                           intercept = 10,
                                           coefs = c(b1, b2),
                                           predictors = list(runif(n, 0, 1),
                                                             rbinom(n, 1, prob = 0.5)),
                                           sigma = 0.7,
                                           interaction_coefs = c(interaction_effect))

#plot, colored by treatment
ggplot(simTwoPredInteraction, aes(x = pred1, y = yobs, color = pred2)) +
  geom_point()

#fit model using rstanarm
fit2 = rstanarm::stan_glm(yobs ~ pred1*pred2, data = simTwoPredInteraction)
coef(fit2)
sigma(fit2)

#for me this yields:
# intercept = 9.8 (true intercept = 10)
# pred1 coefficient = 2.3 (true coefficient = 2)
# pred2 coefficient = 0.08 (true coefficient = 0)
# interaction coefficient = 0.94 (true coefficient = 1)
#sigma = 0.6 (true sigma = 0.7)
#not horrible

#when I vary the sigma (e.g., sigma = 0.1) these become MUCH more accurate, 
# (e.g., a = 9.9, b1 = 2.0, b2 = 0.02, b3 = 0.96) so I'm confident the code works


#####################################
#### PART 2
####################################
#vary sample size
#using just 20% of your data, see how well you can return the parameters you set
n = 20
#for a single continuous predictor:
simOnePred20 = simulateLinRegData(sample_size = n,
                                  intercept = 10,
                                  coefs = c(2),
                                  predictors = list(runif(n, 0, 1)),
                                  sigma = 0.7)

#fit model
fit20 = rstanarm::stan_glm(yobs~pred1, data = simOnePred20)
coef(fit20)
#for me this yields: intercept = 10.3 (true = 10), predictor 1 coeff = 1.97 (true = 2)
#this is still surprisingly accurate for only 20 datapoints

#plot
#plot with fit line
ggplot(simOnePred20, aes(x = pred1, y = yobs)) +
  geom_point() +
  geom_abline(intercept = coef(fit20)[1], slope = coef(fit20)[2], col = 'red')


#now with two predictors and their interaction!
b1 = 2
b2 = 0
interaction_effect = b1/2

simTwoPredInter20 = simulateLinRegData(sample_size = n,
                                           intercept = 10,
                                           coefs = c(b1, b2),
                                           predictors = list(runif(n, 0, 1),
                                                             rbinom(n, 1, prob = 0.5)),
                                           sigma = 0.7,
                                           interaction_coefs = c(interaction_effect))

#plot, colored by treatment
ggplot(simTwoPredInter20, aes(x = pred1, y = yobs, color = pred2)) +
  geom_point()

#fit model using rstanarm
fit_inter20 = rstanarm::stan_glm(yobs ~ pred1*pred2, data = simTwoPredInter20)
coef(fit_inter20)

#for me this yields:
#intercept = 10.3 (true = 10)
#b1 = 1.6 (true = 2)
#b2 = 0.4 (true = 0)
#interaction = -0.3 (true = 1)

#yikes that's bad!! unsurprisingly, you can't accurately fit a complex model
#with a very small sample size. Seems to get worse the higher order the interaction
#e.g., intercept/sigma not bad, slopes mediocre, interaction BAD


#vary sample size from 1:100 and plot parameter estimates

#tragically I have to start at 2 (not 1), because rstanarm does NOT like it when
#I try to fit a model off of a single data point. go figure.
n = c(2:100)

#FIRST, A SINGLE CONTINUOUS PREDICTOR
#initialize a dataframe to hold the parameter values that we generate during loop
single_df = data.frame(matrix(nrow = 99, ncol = 4))
colnames(single_df) = c("sample_size", "sigma", "intercept", "b1")

#simulate data
simOnePred = simulateLinRegData(sample_size = 100,
                                intercept = 10,
                                coefs = c(2),
                                predictors = list(runif(100, 0, 1)),
                                sigma = 0.7)

#loop through 2:100 -- for each loop, take a random draw of simulated data, 
#fit a model to that subset, and save the parameters in the dataframe
for (i in n){
  #random draw of size i
  current_sample = simOnePred %>% sample_n(i)
  
  #fit model
  fit = rstanarm::stan_glm(yobs ~ pred1, data = current_sample)
  
  #record results
  single_df[i-1,1] = i
  single_df[i-1, 2] = sigma(fit)
  single_df[i-1, 3] = coef(fit)[1]
  single_df[i-1, 4] = coef(fit)[2]
}

#got a few divergent transitions & such, presumably for the models with VERY FEW
#datapoints. again, go figure.

#make a grid plot of different parameters
parameters = c("sigma", "intercept", "b1")
true_values = c(0.7, 10, 2)
plot.list.1 = list()

for (i in 1:length(parameters)){
  p = eval(substitute(
    ggplot(single_df, aes(x = sample_size, y = single_df[,i+1], color = "salmon")) +
      geom_point() +
      geom_hline(aes(yintercept = true_values[i])) +
      labs(y = "Parameter Estimate", x = "Sample Size", title = parameters[i]) +
      theme_classic() +
      theme(legend.position="none"),
    list(i = i)))
  plot.list.1[[i]] = p
}

singlepredictor_grid = grid.arrange(plot.list.1[[1]], 
                                    plot.list.1[[2]], 
                                    plot.list.1[[3]],
                                    ncol = 1)
ggsave("figures/part2/singlepredictor_grid.jpg",
       singlepredictor_grid,
       width = 1000, height = 3000,
       units = c("px"))

#To me, it looks like the estimates for the slopes are a little bit worse at 
#small sample sizes than the sigma and intercept estimates (there is a wider 
#spread for b1 at n < 25, despite the fact that its true value (2) is smaller 
#than the true intercept value (10)). 
#Mathematically speaking I'm not sure why this occurs. Google tells me it is 
#because outliers have more of an effect on slope than on intercept. This seems 
#like a reasonable explanation to me, although my intuition here is a little fuzzy.


# NEXT, TWO PREDICTORS AND THEIR INTERACTION!
#make an empty df to store parameters (sample_size, sigma, intercept, b1, b2, b3)
interaction_df = data.frame(matrix(nrow = 99, ncol = 6))
colnames(interaction_df) = c("sample_size", "sigma", "intercept", "b1", "b2", "b3")

#set parameters for simulation
b1 = 2
b2 = 0
b3 = b1/2

#simulate data
simTwoPredInter = simulateLinRegData(sample_size = 100,
                                       intercept = 10,
                                       coefs = c(b1, b2),
                                       predictors = list(runif(100, 0, 1),
                                                         rbinom(100, 1, prob = 0.5)),
                                       sigma = 0.7,
                                       interaction_coefs = c(b3))

#loop through 2:n, taking a random draw of size i from simulated data, fitting
#the model and recording parameter values for each iteration
for (i in n){
  
  #for small n, rstanarm will throw an error if predictor 2 (binary)
  #only take a single value. Therefore, draw half the values from each group
  #BUT, only do this up until n = 50. At a certain point, the probability of 
  #only getting 1 OR 0 is infinitismally small. But we don't want to do this
  #for very large n (e.g., n = 100) because there likely isn't a perfect 50/50
  #split in the original data (and it will throw an error)
  if( i < 50){
  current_sample = simTwoPredInter %>%
    group_by(pred2) %>%
    sample_n( i%/% 2 + if_else(i%%2 > 0 & cur_group_id() == 1, 1, 0))
  }
  else{
    current_sample = simTwoPredInter %>% sample_n(i)
  }
  
  #fit model
  fit = rstanarm::stan_glm(yobs ~ pred1*pred2, data = current_sample)
  
  #record results
  interaction_df[i-1,1] = i
  interaction_df[i-1, 2] = sigma(fit)
  interaction_df[i-1, 3] = coef(fit)[1]
  interaction_df[i-1, 4] = coef(fit)[2]
  interaction_df[i-1, 5] = coef(fit)[3]
  interaction_df[i-1, 6] = coef(fit)[4]
}

#more divergent transitions on this one--a couple of the models didn't actually
#converge due to small sample size

#plot results
parameters = c("sigma", "intercept", "b1", "b2", "b3")
true_values = c(0.7, 10, 2, 0, 1)
plot.list.2 = list()

for (i in 1:length(parameters)){
  p = eval(substitute(
    ggplot(interaction_df, aes(x = sample_size, y = interaction_df[,i+1], color = "salmon")) +
      geom_point() +
      geom_hline(aes(yintercept = true_values[i])) +
      labs(y = "Parameter Estimate", x = "Sample Size", title = parameters[i]) +
      theme_classic() +
      theme(legend.position="none"),
    list(i = i)))
  plot.list.2[[i]] = p
}

interactingpredictors_grid = grid.arrange(plot.list.2[[1]], 
                                    plot.list.2[[2]], 
                                    plot.list.2[[3]],
                                    plot.list.2[[4]],
                                    plot.list.2[[5]],
                                    ncol = 2)
ggsave("figures/part2/interactingpredictors_grid.jpg",
       interactingpredictors_grid,
       width = 2000, height = 3000,
       units = c("px"))

#well...sigma is certainly more accurately estimated than the other parameters!
#From these plots you can see that the full initial dataset (n = 100) doesn't
#even lead to very accurate estimates for the interaction coefficient or the 
#predictor 2 coefficient. this could be in part due to an "unlucky" simulation,
#so the monte carlo in the next step may improve it. But this follows a similar trend
#from the model with only a single intercept: the intercept and sigma parameters
#are more accurately fit than the slope parameters! Of the slopes, I would say
#that the interaction estimate is the worst estimate (largest variation in estimates
#and it "converges" (settles?) on a value near n=100 that is quite far from the
#actual simulating value)


#Improve sampling with monte carlo simulation of data
#Alrighty let's see if fitting 990 models makes my macbook air explode!

#FIRST, A SINGLE CONTINUOUS PREDICTOR
#set sample size range
n = c(2:100)

#initialize df to hold parameter estimates
single_df_mc = data.frame(matrix(nrow = 990, ncol = 4))
colnames(single_df_mc) = c("sample_size", "sigma", "intercept", "b1")

#set a count variable (to tell us which iteration we're on, and add parameter
#estimates to the proper row in the dataframe)
count = 1

#nested for loop. outer loop: same as above, loops over sample sizes from 2-100
#inner loop repeats each sample size value 10 times
#this time the data simulation occurs inside the loop
for (i in n){
  for (j in 1:10){
  #simulate data
  simOnePred = simulateLinRegData(sample_size = i,
                                  intercept = 10,
                                  coefs = c(2),
                                  predictors = list(runif(i, 0, 1)),
                                  sigma = 0.7)
  
  #fit model
  fit = rstanarm::stan_glm(yobs ~ pred1, data = simOnePred)
  
  #record results
  single_df_mc[count,1] = i
  single_df_mc[count, 2] = sigma(fit)
  single_df_mc[count, 3] = coef(fit)[1]
  single_df_mc[count, 4] = coef(fit)[2]
  count = count + 1
  
  }
}

#plot results
#make a grid plot of different parameters
parameters = c("sigma", "intercept", "b1")
true_values = c(0.7, 10, 2)
plot.list.3 = list()

for (i in 1:length(parameters)){
  p = eval(substitute(
    ggplot(single_df_mc, aes(x = sample_size, y = single_df_mc[,i+1], color = "salmon")) +
      geom_point() +
      geom_hline(aes(yintercept = true_values[i])) +
      labs(y = "Parameter Estimate", x = "Sample Size", title = parameters[i]) +
      scale_y_continuous(limits=c(true_values[i]-5, true_values[i] + 5)) +
      theme_classic() +
      theme(legend.position="none"),
    list(i = i)))
  plot.list.3[[i]] = p
}

singlepredictor_mc_grid = grid.arrange(plot.list.3[[1]], 
                                    plot.list.3[[2]], 
                                    plot.list.3[[3]],
                                    ncol = 1)
ggsave("figures/part2/singlepredictor_mc_grid.jpg",
       singlepredictor_mc_grid,
       width = 1000, height = 3000,
       units = c("px"))

#In order of fit quality, I would say sigma > intercept > slope, especially
#at small sample sizes. again keeping with previous patterns



# NEXT, TWO PREDICTORS AND THEIR INTERACTION!

#initialize df to store parameter estimates (sample_size, sigma, intercept, b1, b2, b3)
interaction_df_mc = data.frame(matrix(nrow = 990, ncol = 6))
colnames(interaction_df_mc) = c("sample_size", "sigma", "intercept", "b1", "b2", "b3")

#set parameters for simulation
b1 = 2
b2 = 0
b3 = b1/2

#set count variable to track iteration
count = 1

#effectively the same for loop as last part, but for more complex simulation/model
for (ii in n){
  for(jj in 1:10){
  #simulate data
  simTwoPredInter = simulateLinRegData(sample_size = ii,
                                       intercept = 10,
                                       coefs = c(b1, b2),
                                       predictors = list(runif(ii, 0, 1),
                                                         rbinom(ii, 1, prob = 0.5)),
                                       sigma = 0.7,
                                       interaction_coefs = c(b3))
  
  #check for unique values of pred 2, and coerce to different values if constant
  if(length(unique(simTwoPredInter$pred2)) == 1){
    simTwoPredInter$pred2[1] = 1
    simTwoPredInter$pred2[2] = 0
  }

  #fit model
  fit = rstanarm::stan_glm(yobs ~ pred1*pred2, data = simTwoPredInter)
  
  #record results
  interaction_df_mc[count,1] = ii
  interaction_df_mc[count, 2] = sigma(fit)
  interaction_df_mc[count, 3] = coef(fit)[1]
  interaction_df_mc[count, 4] = coef(fit)[2]
  interaction_df_mc[count, 5] = coef(fit)[3]
  interaction_df_mc[count, 6] = coef(fit)[4]
  count = count + 1
  }
}

#plot results
parameters = c("sigma", "intercept", "b1", "b2", "b3")
true_values = c(0.7, 10, 2, 0, 1)
plot.list.4 = list()

for (i in 1:length(parameters)){
  p = eval(substitute(
    ggplot(interaction_df_mc, aes(x = sample_size, y = interaction_df_mc[,i+1], color = "salmon")) +
      geom_point() +
      geom_hline(aes(yintercept = true_values[i])) +
      labs(y = "Parameter Estimate", x = "Sample Size", title = parameters[i]) +
      scale_y_continuous(limits=c(true_values[i]-5, true_values[i]+5)) +
      theme_classic() +
      theme(legend.position="none"),
    list(i = i)))
  plot.list.4[[i]] = p
}

interactingpredictors_mc_grid = grid.arrange(plot.list.4[[1]], 
                                          plot.list.4[[2]], 
                                          plot.list.4[[3]],
                                          plot.list.4[[4]],
                                          plot.list.4[[5]],
                                          ncol = 2)
ggsave("figures/part2/interactingpredictors_mc_grid.jpg",
       interactingpredictors_mc_grid,
       width = 2000, height = 3000,
       units = c("px"))


#####################################
#### PART 3
####################################
#everything exactly the same but we're changing sigma to 1.05 (50% increase)

#FIRST, A SINGLE CONTINUOUS PREDICTOR
#set sample size range
n = c(2:100)

#initialize df to hold parameter estimates
single_df_mc_higherror = data.frame(matrix(nrow = 990, ncol = 4))
colnames(single_df_mc_higherror) = c("sample_size", "sigma", "intercept", "b1")

#set a count variable (to tell us which iteration we're on, and add parameter
#estimates to the proper row in the dataframe)
count = 1

#nested for loop. outer loop: same as above, loops over sample sizes from 2-100
#inner loop repeats each sample size value 10 times
#this time the data simulation occurs inside the loop
for (i in n){
  for (j in 1:10){
    #simulate data
    simOnePred = simulateLinRegData(sample_size = i,
                                    intercept = 10,
                                    coefs = c(2),
                                    predictors = list(runif(i, 0, 1)),
                                    sigma = 1.05)
    
    #fit model
    fit = rstanarm::stan_glm(yobs ~ pred1, data = simOnePred)
    
    #record results
    single_df_mc_higherror[count,1] = i
    single_df_mc_higherror[count, 2] = sigma(fit)
    single_df_mc_higherror[count, 3] = coef(fit)[1]
    single_df_mc_higherror[count, 4] = coef(fit)[2]
    count = count + 1
    
  }
}

#plot results
#make a grid plot of different parameters
parameters = c("sigma", "intercept", "b1")
true_values = c(1.05, 10, 2)
plot.list.5 = list()

for (i in 1:length(parameters)){
  p = eval(substitute(
    ggplot(single_df_mc_higherror, aes(x = sample_size, y = single_df_mc_higherror[,i+1], color = "salmon")) +
      geom_point() +
      geom_hline(aes(yintercept = true_values[i])) +
      labs(y = "Parameter Estimate", x = "Sample Size", title = parameters[i]) +
      scale_y_continuous(limits = c(true_values[i]-5, true_values[i]+5)) +
      theme_classic() +
      theme(legend.position="none"),
    list(i = i)))
  plot.list.5[[i]] = p
}

singlepredictor_mc_higherror_grid = grid.arrange(plot.list.5[[1]], 
                                       plot.list.5[[2]], 
                                       plot.list.5[[3]],
                                       ncol = 1)
ggsave("figures/part3/singlepredictor_mc_higherror_grid.jpg",
       singlepredictor_mc_higherror_grid,
       width = 1000, height = 3000,
       units = c("px"))

# NEXT, TWO PREDICTORS AND THEIR INTERACTION!

#initialize df to store parameter estimates (sample_size, sigma, intercept, b1, b2, b3)
interaction_df_mc_higherror = data.frame(matrix(nrow = 990, ncol = 6))
colnames(interaction_df_mc_higherror) = c("sample_size", "sigma", "intercept", "b1", "b2", "b3")

#set parameters for simulation
b1 = 2
b2 = 0
b3 = b1/2

#set count variable to track iteration
count = 1

#effectively the same for loop as last part, but for more complex simulation/model
for (ii in n){
  for(jj in 1:10){
    #simulate data
    simTwoPredInter = simulateLinRegData(sample_size = ii,
                                         intercept = 10,
                                         coefs = c(b1, b2),
                                         predictors = list(runif(ii, 0, 1),
                                                           rbinom(ii, 1, prob = 0.5)),
                                         sigma = 1.05,
                                         interaction_coefs = c(b3))
    
    #check for unique values of pred 2, and coerce to different values if constant
    if(length(unique(simTwoPredInter$pred2)) == 1){
      simTwoPredInter$pred2[1] = 1
      simTwoPredInter$pred2[2] = 0
    }
    
    #fit model
    fit = rstanarm::stan_glm(yobs ~ pred1*pred2, data = simTwoPredInter)
    
    #record results
    interaction_df_mc_higherror[count,1] = ii
    interaction_df_mc_higherror[count, 2] = sigma(fit)
    interaction_df_mc_higherror[count, 3] = coef(fit)[1]
    interaction_df_mc_higherror[count, 4] = coef(fit)[2]
    interaction_df_mc_higherror[count, 5] = coef(fit)[3]
    interaction_df_mc_higherror[count, 6] = coef(fit)[4]
    count = count + 1
  }
}

#plot results
parameters = c("sigma", "intercept", "b1", "b2", "b3")
true_values = c(1.05, 10, 2, 0, 1)
plot.list.6 = list()

for (i in 1:length(parameters)){
  p = eval(substitute(
    ggplot(interaction_df_mc_higherror, aes(x = sample_size, y = interaction_df_mc_higherror[,i+1], color = "salmon")) +
      geom_point() +
      geom_hline(aes(yintercept = true_values[i])) +
      labs(y = "Parameter Estimate", x = "Sample Size", title = parameters[i]) +
      scale_y_continuous(limits=c(true_values[i]-5, true_values[i]+5)) +
      theme_classic() +
      theme(legend.position="none"),
    list(i = i)))
  plot.list.6[[i]] = p
}

interactingpredictors_mc_higherror_grid = grid.arrange(plot.list.6[[1]], 
                                             plot.list.6[[2]], 
                                             plot.list.6[[3]],
                                             plot.list.6[[4]],
                                             plot.list.6[[5]],
                                             ncol = 2)
ggsave("figures/part3/interactingpredictors_mc_higherror_grid.jpg",
       interactingpredictors_mc_higherror_grid,
       width = 2000, height = 3000,
       units = c("px"))

#####################################
#### CHALLENGE 1
####################################
#I'm going to stop using generalizable functions for simulation here because I'm lazy

#check out function simulateChallenge1 for simulation specification

simDataChallenge1 = simulateChallenge1(sample_size = 1000,
                                       intercept = 10,
                                       sigma = 0.7,
                                       b1 = 2,
                                       btreatment1 = 1,
                                       btreatment2 = 2,
                                       interaction_treatment1 = 0.5,
                                       interaction_treatment2 = 1.5)

#fit models
fit.dummy = rstanarm::stan_glm(yobs ~ pred1*treatment1 + pred1*treatment2, data = simDataChallenge1)
summary(fit.dummy)
#the predictions for n = 100 are not fantastic, but I've checked the code with n = 1000
#and it's quite close!

#going to go straight to the monte carlo simulation bit!
params = monteCarloChallenge1(number_repeats = 10,
                              sample_size = 100,
                              intercept = 10,
                              sigma = 0.7,
                              b1 = 2,
                              btreatment1 = 1,
                              btreatment2 = 2,
                              interaction_treatment1 = 0.5,
                              interaction_treatment2 = 1.5
  
)

#plotting results
plots = colnames(params)[-1]
true_values = c(0.7, 10, 2, 1, 2, 0.5, 1.5)
plot.list = list()

for (i in 1:length(plots)){
  p = eval(substitute(
    ggplot(params, aes(x = sample_size, y = params[,i+1], color = "salmon")) +
    geom_point() +
    geom_hline(aes(yintercept = true_values[i])) +
    labs(y = "Parameter Estimate", x = "Sample Size", title = plots[i]) +
    scale_y_continuous(limits=c(true_values[i]-5, true_values[i]+5)) +
    theme_classic() +
    theme(legend.position="none"),
    list(i = i)))
  print(i)
  print(p)
  plot.list[[i]] = p
}

challenge1grid = grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]],
                               plot.list[[4]], plot.list[[5]], plot.list[[6]],
                               plot.list[[7]], ncol = 2)
ggsave("figures/challenges/challenge1grid.jpg",
       challenge1grid,
       width = 1500, height = 3000,
       units = c("px"))



#####################################
#### CHALLENGE 2 & 3
####################################
#two binary treatments (e.g., control/treatment + male/female)
#two and three way interactions
#going to do these challenges simultaneously because I don't want to run the code twice!

#simulate data
threePredSim = simulateChallenge2(sample_size = 100,
                                 intercept = 10,
                                 sigma = 0.7,
                                 beta = c(2, 3, 4, 1, 1.5, 2, 0.5))

#run model
fit3 = rstanarm::stan_glm(yobs ~ pred1*pred2*pred3, data = threePredSim)
summary(fit3)

#monte carlo
paramschallenge2 = monteCarloChallenge2(number_repeats = 10,
                              sample_size = 100,
                              intercept = 10,
                              sigma = 0.7,
                              beta = c(2, 3, 4, 1, 1.5, 2, 0.5),
                              confidence_interval = 0.5
                              
)

#plot results
sigma = ggplot(paramschallenge2, aes(x=sample_size, y=sigma)) + 
  geom_point() +
  geom_errorbar(aes(ymin=sigmalower, ymax=sigmaupper)) +
  geom_hline(yintercept = 0.7) +
  scale_y_continuous(limits=c(0.7-20, 0.7+20)) +
  labs(x = "Sample Size", y = "Sigma Estimate")
intercept = ggplot(paramschallenge2, aes(x=sample_size, y=intercept)) + 
  geom_point() +
  geom_errorbar(aes(ymin=interceptlower, ymax=interceptupper)) +
  geom_hline(yintercept = 10) +
  scale_y_continuous(limits=c(10-20, 10+20)) +
  labs(x = "Sample Size", y = "Intercept Estimate")
b1 = ggplot(paramschallenge2, aes(x=sample_size, y=b1)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b1lower, ymax=b1upper)) +
  geom_hline(yintercept = 2) +
  scale_y_continuous(limits=c(2-20, 2+20)) +
  labs(x = "Sample Size", y = "B1 Estimate")
b2 = ggplot(paramschallenge2, aes(x=sample_size, y=b2)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b2lower, ymax=b2upper)) +
  geom_hline(yintercept = 3) +
  scale_y_continuous(limits=c(3-20, 3+20)) +
  labs(x = "Sample Size", y = "B2 Estimate")
b3 = ggplot(paramschallenge2, aes(x=sample_size, y=b3)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b3lower, ymax=b3upper)) +
  geom_hline(yintercept = 4) +
  scale_y_continuous(limits=c(4-20, 4+20)) +
  labs(x = "Sample Size", y = "B3 Estimate")
b4 = ggplot(paramschallenge2, aes(x=sample_size, y=b4)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b4lower, ymax=b4upper)) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(limits=c(1-20, 1+20)) +
  labs(x = "Sample Size", y = "B4 (Pred 1*Pred 2) Estimate")
b5 = ggplot(paramschallenge2, aes(x=sample_size, y=b5)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b5lower, ymax=b5upper)) +
  geom_hline(yintercept = 1.5) +
  scale_y_continuous(limits=c(1.5-20, 1.5+20)) +
  labs(x = "Sample Size", y = "B5 (Pred 1*Pred 3) Estimate")
b6 = ggplot(paramschallenge2, aes(x=sample_size, y=b6)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b6lower, ymax=b6upper)) +
  geom_hline(yintercept = 2) +
  scale_y_continuous(limits=c(2-20, 2+20)) +
  labs(x = "Sample Size", y = "B6 (Pred 2*Pred 3) Estimate")
b7 = ggplot(paramschallenge2, aes(x=sample_size, y=b7)) + 
  geom_point() +
  geom_errorbar(aes(ymin=b7lower, ymax=b7upper)) +
  geom_hline(yintercept = 0.5) +
  scale_y_continuous(limits=c(0.5-20, 0.5+20)) +
  labs(x = "Sample Size", y = "B7 (Pred1*Pred2*Pred3) Estimate")

threewayinteraction_grid = grid.arrange(sigma,
                                        intercept,
                                        b1,
                                        b2,
                                        b3,
                                        b4,
                                        b5,
                                        b6,
                                        b7,
                                        ncol = 3)
ggsave("figures/challenges/threewayinteraction_grid.jpg",
       threewayinteraction_grid,
       width = 3000, height = 3000,
       units = c("px"))




#wowza look at the yscale on that 3 way interaction estimate (some really big errors!!)