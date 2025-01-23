#Abi Brown Homework 1 for Jan 2025 Bayesian modelling

#packages and options
library(rstanarm)
options(stringsasfactors = FALSE)
# Part 1 ------------------------------------------------------------------

#building a linear model and simulating data
n <- 100
sig <- 4
b <- 4
a <- 10
xdata <- runif(n, min = 0, max = 12)
ypred <- a+b*xdata
yobs <- ypred + rnorm(n, 0, sig)
model1 <- stan_glm(yobs~xdata)
print(summary(model1))
#model is able to return good estimates of parameters
#on the first run mean estimates were sig = 3.4, b = 3.9, a = 10.4

#adding binary covariate with an interaction
covdata <- rbinom(n, size = 1, prob = 0.5)  
#effect of treatment
b1 <- 1
#interaction with continuous variable
b2 <- 2
#new model with interaction
ypred2 <- a+b*xdata+b1*covdata+b2*xdata*covdata
yobs2 <- ypred2 +rnorm(n, 0, sig)
model2 <- stan_glm(yobs2~xdata+covdata+xdata*covdata)
print(summary(model2))
#model was able to return most parameters, although estimates are worse than without covariate
#estimate of the covariate is quite bad, I am not sure if this is because of the model (too small of an effect size)
##or because I made an error in simulating the data
#first run estimates mean estimates:a=12.1, b=3.6, b1=-1.7, b2=2.6, sig=4.1) 

# Part 2 ------------------------------------------------------------------

#returning parameters with only 20% of the data
n <- 20
xdata.s <- sample(xdata,n)
covdata.s <- sample(covdata,n)
ypred3 <- a+b*xdata.s+b1*covdata.s+b2*xdata.s*covdata.s
yobs3 <- ypred3+rnorm(n,0,sig)
model3 <- stan_glm(yobs3~xdata.s+covdata.s+xdata.s*covdata.s)
print(summary(model3))
#estimates of b (first variable slope), interaction, and error were still close, but binary covariate slope estimate got much worse

#estimated parameters for 1:n of of our data
for(n in c(5,10,20,30,40,50,60,70,80,90,100)){
  xdata.s <- sample(xdata,n)
  covdata.s <- sample(covdata,n)
  yobs4 <- a+b*xdata.s+b1*covdata.s+b2*xdata.s*covdata.s+rnorm(n,0,sig)
  model <- stan_glm(yobs4~xdata.s+covdata.s+xdata.s*covdata.s)
  name <- paste0("model", n)
  assign(name, as.data.frame(t(as.data.frame(summary(model)))))
}  
#sample size of 5 was too low, model run received many errors (sample size and divergent transitions),
#plots start at sample size 10
#data frame for parameter estimates at each sample size
#adding sample size column
model10$n <- 10
model20$n <- 20
model30$n <- 30
model40$n <- 40
model50$n <- 50
model60$n <- 60
model70$n <- 70
model80$n <- 80
model90$n <- 90
model100$n <- 100
param.ests <- rbind.data.frame(model10, model20, model30, model40, model50, model60, model70, model80, model90, model100)
#removing non-means since we wont use them for plotting
param.ests <- param.ests[-(2:8),]
param.ests <- param.ests[-(3:9),]
param.ests <- param.ests[-(4:10),]
param.ests <- param.ests[-(5:11),]
param.ests <- param.ests[-(6:12),]
param.ests <- param.ests[-(7:13),]
param.ests <- param.ests[-(8:14),]
param.ests <- param.ests[-(9:15),]
param.ests <- param.ests[-(10:16),]
param.ests <- param.ests[-(11:17),]
#changing column names so I can plot
colnames(param.ests) <- c("intercept", "b", "b1", "b2", "sigma", "mean_PPD", "log_posterior", "n")
#I was having trouble saving plot it code, so I have downloaded them manually in RStudio
plot(param.ests$n, param.ests$intercept, xlab = "n", ylab = "Intercept estimate", main = "Intercept - true value  = 10")
plot(param.ests$n, param.ests$b, xlab = "n", ylab = "continuous var effect estimate", main = "b - true value = 4")
plot(param.ests$n, param.ests$b1, xlab = "n", ylab = "binary var effect estimate", main = "b1 (binary covariate) - true value = 1")
plot(param.ests$n, param.ests$b2, xlab = "n", ylab = "interaction effect estimate", main = "b2 (interaction) - true value = 2")
plot(param.ests$n, param.ests$sig, xlab = "n", ylab = "error estimate", main = "sigma - true value = 2")

#For most parameters, estimates tend to get closer to the true value as sample size increases
#we are taking more samples, so there is less chance of obtaining a biased sample and we are getting closer to representing the true distribution
#estimates of the binary variable effect error bounce around a lot, even as sample numbers increase
#I think this might have something to do with how I simulated the binary covariate data, especially at low sample sizes
##where there is a high chance that there are few or no samples from either of the treatment categories (data is all 0s or all 1s)


#Improving sampling 
#I know that this code does not work, I would really like to spend some time going over this in class
#for(n in c(5,10,20,35,50,75,100)){
# for(i in 1:10) {
# xdata[i] <- runif(n, min = 0, max = 12)
#  covdata[i] <- rbinom(n, size =1, prob = 0.5)
#  xdata.s[i] <- sample(xdata[i],n)
# covdata.s[i] <- sample(covdata[i],n)
# yobs4[i] <- a+b*xdata.s[i]+b1*covdata.s[i]+b2*xdata.s[i]*covdata.s[i]+rnorm(n,0,sig)
#  model <- stan_glm(yobs4[i]~xdata.s[i]+covdata.s[i]+xdata.s[i]*covdata.s[i])
#  name <- paste0("model10draw", i, n)
#  assign(name, model)} }
#if I understand this correctly, we want 10 runs of the stan_glm model for each sample size, and for each of these runs to 
##draw from a newly generated simulated dataset
#I am unsure of how we would aggregate/represent the set of model estimates for each sample size after that
#it would be great to go over this in class - more generally about variability between model runs and how to deal with it

#increasing the error term by 50%
sig <- 6
for(n in c(5,10,20,30,40,50,60,70,80,90,100)){
  xdata.s <- sample(xdata,n)
  covdata.s <- sample(covdata,n)
  yobs4 <- a+b*xdata.s+b1*covdata.s+b2*xdata.s*covdata.s +rnorm(n,0,sig)
  model <- stan_glm(yobs4~xdata.s+covdata.s+xdata.s*covdata.s)
  name <- paste0("model", n)
  assign(name, model)
}  
#data frame for parameter estimates at each sample size
#adding sample size column
model10$n <- 10
model20$n <- 20
model30$n <- 30
model40$n <- 40
model50$n <- 50
model60$n <- 60
model70$n <- 70
model80$n <- 80
model90$n <- 90
model100$n <- 100
param.ests2 <- rbind.data.frame(model10, model20, model30, model40, model50, model60, model70, model80, model90, model100)
#removing non-means since we wont use them for plotting
param.ests2<- param.ests2[-(2:8),]
param.ests2 <- param.ests2[-(3:9),]
param.ests2 <- param.ests2[-(4:10),]
param.ests2 <- param.ests2[-(5:11),]
param.ests2 <- param.ests2[-(6:12),]
param.ests2 <- param.ests2[-(7:13),]
param.ests2 <- param.ests2[-(8:14),]
param.ests2 <- param.ests2[-(9:15),]
param.ests2 <- param.ests2[-(10:16),]
param.ests2 <- param.ests2[-(11:17),]
#changing column names so I can plot
colnames(param.ests2) <- c("intercept", "b", "b1", "b2", "sigma", "mean_PPD", "log_posterior", "n")
#plots are downloaded in figures folder with filenames with hisig
plot(param.ests2$n, param.ests2$intercept, xlab = "n", ylab = "Intercept estimate", main = "Intercept - true value  = 10")
plot(param.ests2$n, param.ests2$b, xlab = "n", ylab = "continuous var effect estimate", main = "b - true value = 4")
plot(param.ests2$n, param.ests2$b1, xlab = "n", ylab = "binary var effect estimate", main = "b1 (binary covariate) - true value = 1")
plot(param.ests2$n, param.ests2$b2, xlab = "n", ylab = "interaction effect estimate", main = "b2 (interaction) - true value = 2")
plot(param.ests2$n, param.ests2$sig, xlab = "n", ylab = "error estimate", main = "sigma - true value = 3")
#error estimates became a lot less reliable with the increase
#estmates of other parameters (like b) also got a lot less close to true values even at 100 samples
#this makes sense as increasing the variability in the data makes sampling the same amount from the posterior less likely to capture "true values"

#Some of my estimates at low sample sizes are likely not reliable as models had divergent transitions
#I took a brief look at the documentation for this error but it would be great to go over what this means in class, as well as
#what model non-convergence means