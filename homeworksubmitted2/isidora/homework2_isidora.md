---
title: "Homework2"
author: "Isidora Silva-Valderrama"
date: "2025-01-23"
output: html_document
---

# Homework 2

```{r setup, include=FALSE}
#setwd
setwd("C:/Users/ibsil/Dropbox/ISIDORA (por Isidora)/UBC/FRST507C - Hierarchical mode l& Bayesian stats")

options(mc.cores = parallel::detectCores())
rm(list = ls()) #cleans your R space
options(stringsAsFactor = FALSE) #to stop R from automatically determining that your data is a factor

```

from Lizzie`s dataset: 
analyses/input/carnivores
```{r}
carn <- read.csv("./carnivorebodymass.csv")
teeth <- read.csv("./carnivoreteeth.csv")

#unite both datasets
library(dplyr)
carn <- carn %>%
  rename(Species = species)

data <- merge(carn,teeth, by = c("Species", "sex", "Family"))
data <- data[,c(1:5,8:16)]
```


## Come up with a model

### Exploring the data
```{r}
head(data)
```

Metadata: 
- mass.midpoint = represent median species values derived from multiple sources
- PM4 = premolar / length of the upper carnassial 
- CsupL = canine / length of the lower carnassial? 

I would like to fit a model to predict the mass of a carnivore based on their teeth length. The data was obtained from a museum archives, so I would think that it could have relevance for archeological findings, for example.
In that case, I would think that the length of the canine is positively correlated with the body mass of the carnivore. Bigger animals would have bigger teeth.

```{r}
length(unique(data$Species))

table(data$sex)

#visualization
library(ggplot2)
ggplot(data, aes(x=CsupL))+  #my x
  geom_histogram()

ggplot(data, aes(x=CsupL, y=mass.midpoint))+ # x and y relation
  geom_point()+
  theme_classic()

ggplot(data, aes(x=log(CsupL), y=log.midpoint.mass))+ # x and y relation
  geom_point()+
  theme_classic()
```

So there is definitely a trend between the mass of the carnivore and the teeth length, and we have different groups: species, genus or family. Just for simplicity (due to the time constrains for the homework), I will make a model that can predict the mass of a carnivore from the length of its upper canines. 

For a simple linear model:

y_simulated = alpha + beta1*x1 + error

body_mass = alpha + beta1*canine_length + error

The distribution of the canine data looks different than a normal gaussian curve, as it is skewed to the left. I will then simulate the data following a lognormal distribution. 

## 2. Simulate data
```{r}
n <- 1000

CsupL_sim <- rlnorm(n, meanlog = 1.4, sdlog = .7)  # Left-skewed (x1)
hist(CsupL_sim)

#my parameters
alpha <- 1
beta1 <- 1.3
sigma <- 0.3

error <- rnorm(n,0,sigma)

log.midpoint.mass_sim <- alpha + beta1*log(CsupL_sim) + error
plot(log.midpoint.mass_sim ~ log(CsupL_sim), xlim = c(0,4), ylim = c(2,6))
```

I can also add the sex as a binary predictor. The body mass of the carnivore can be different for females and males. I would say that females are usually smaller than males. 

body_mass = alpha + beta1 x canine_length + beta2 x sex + error

```{r}
sex_sim <- rbinom(n, size = 1, prob = c(0.5,0.5)) #prob: half male, half female (x2)
beta2 <- 0.5

log.midpoint.mass_sim <- alpha + beta1*log(CsupL_sim) + beta2*sex_sim + error
```

Fit the model with simulated data: 

```{r}
library(rstanarm)
mod_sim <- stan_glm(log.midpoint.mass_sim ~ log(CsupL_sim) + sex_sim,
                      seed= 1234)
summary(mod_sim)

```

From the summary above, I could retrieve all my parameters. So, the model does not have mathematical incongruities or mistakes. 

## 3. Develop priors

Because rstanarm uses default prior distributions that are weakly informative, I don't need to specify them before but it is always useful to know what was used and to make them explicit. 
```{r}
prior_summary(object = mod_sim)
```
## 4. Run your model on empirical data

To practice for future models written in Stan, I'll transform my data into zeros and ones instead of factors: 
```{r}
#Ensure sex a factor
data$sex <- as.factor(data$sex)

#Convert sex to 0/1 encoding
data$sex_binary <- as.numeric(data$sex) - 1
```

Now, female = 0, male = 1. 

```{r}
mod1 <- stan_glm(log.midpoint.mass~log(CsupL) + sex_binary, data=data)
summary(mod1)
```

## 5. Retrodictive checks

```{r}
pairs(mod1)
```
```{r}
#extract posterior draws from the model
posterior <- as.array(mod1) #using as.array instead of as.matrix to keep the markov chains separated.
dimnames(posterior)

#Check the chains
bayesplot::mcmc_trace(posterior) + 
  xlab("Post-warmup iteration")
```

**Posterior uncertainty intervals**
(sometimes called "credible intervals")
```{r}
#Central posterior uncertainty intervals
bayesplot::mcmc_intervals(posterior)
bayesplot::mcmc_areas(posterior)
```

The model does not have any warnings. The chains look good and there are no divergences either. Interestingly, the sex appears to not have great influence on the mass on the log. mass of the carnivore. This was different that what I simulated.

```{r}
pp_check(mod1, ndraws = 1000) + 
  coord_cartesian(xlim = c(1.5, 6)) + ggtitle("Posterior prediction")
```

However, I get a posterior that doen't fit well with the observed data. The multimodality of the curve makes sense because I am not including all the variables that could be affecting the body mass. For example, grouping the effects by genus or species might explain the multimodal posterior that we get. I should, therefore, include it in the model.

## 6. Go back to the model 
Adding Species as a grouping factor:

body_mass <- alpha[species] + beta1[species] x canine_lenght + beta2[species] x sex + error

```{r}
nspecies <- 170 #number of species

Species_sim <- sample(1:nspecies, n, replace=TRUE) 

#my parameters
mu_alpha <- 1
sigma_alpha <- .3

mu_beta1 <- 1.3
sigma_beta1 <- .2

mu_beta2 <- 0.5
sigma_beta2 <- .4

#hyper priors
alpha <- rnorm(nspecies, mu_alpha, sigma_alpha) #intercepts
beta1 <- rnorm(nspecies, mu_beta1, sigma_beta1) #slopes
beta2 <- rnorm(nspecies, mu_beta2, sigma_beta2) #slopes

log.midpoint.mass_sim <- alpha[Species_sim ]  + beta1[Species_sim ]*log(CsupL_sim) + beta2[Species_sim ]*sex_sim
plot(log.midpoint.mass_sim~log(CsupL_sim), xlim = c(0,4), ylim = c(2,6))
```

Fit the model with simulated data: 

```{r}
library(rstanarm)
mod2_sim <- stan_lmer(log.midpoint.mass_sim ~ log(CsupL_sim) + sex_sim + (1|Species_sim),
                      warmup=3000,
                      iter=5000,
                      seed= 1234)
summary(mod2_sim)

```

No warnings. I can retrieve the parameters without problems. Now, with the observed data: the model with Species as a grouping effect with only the intercept
```{r}
mod2 <- stan_lmer(log.midpoint.mass~log(CsupL) + sex_binary + (1|Species), data=data,
                   warmup=3000,
                      iter=6000,
                      seed= 1234,
                  control = list(max_treedepth = 20))
summary(mod2)

```

I still get the warning of too low ESS. I think I might be using too little informative priors. For a next steps, I would try to constrain the priors and see if I still have the problem. 

