# Victor, 18 Jan

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# setup
library(rstan)
library(ggplot2)
wd <- '/home/victor/projects/ubc_related/courses/lizzie_course/bayes2025homework/homeworksubmitted/victor/hmw1'
set.seed(2558)

### question 1.1 ###

# 0. study design
## let's say we have a certain measure of the fitness of a population based on an environmental variable (GDD?)
## in our model, the response of fitness (c[0,1]) to GDD follows a symmetric "bell curve" shape

gaussian <- function(x,a,b,c){
  return(0.7*exp(-(a*x-b)^2/(2*c^2)))
}

# 1. build the model, look at priors
model_file <- file.path(wd, 'model1.stan')
{
  priors <- list()
  priors$a <- rnorm(10000, 1, 0.5)
  priors$b <- rnorm(10000, 2000, 250)
  priors$c <- rnorm(10000, 500, 130)
  priors$sigma <- rnorm(10000, 0, 0.2)
  par(mfrow=c(1, 4))
  hist(priors$a, xlab="", prob=TRUE, 
       xlim=c(-1,2.5), 
       border = FALSE)
  hist(priors$b, xlab="", prob=TRUE, 
       xlim=c(1e3,3e3), 
       border = FALSE)
  hist(priors$c, xlab="", prob=TRUE, 
       xlim=c(0e3,1e3), 
       border = FALSE)
  hist(priors$sigma, xlab="", prob=TRUE, 
       xlim=c(-1,1), 
       border = FALSE)
}

# 2. simulate data, check the model
{
  a <- 1
  n <- 100
  b <- 1700 # optimal GDD for the population
  c <- 600
  sigma <- 0.05
  x <- runif(n, 0, 4000) # simulate some annual GDD, typically <4000?
  yhat <- 0.1 + gaussian(x,a,b,c)
  y <- yhat + rnorm(n, mean = 0, sd = sigma)
  par(mfrow=c(1, 1))
  plot(y ~ x, frame.plot = FALSE)
}
fit  <- stan(model_file, data=list(N=n, y= y, gdd = x), iter=2000, chains=4, cores = 4)
posteriors <- extract(fit, pars = c("a","b", "c", 'sigma'))
{
  par(mfrow=c(1, 4))
  hist(posteriors$a, xlab="", prob=TRUE, 
       xlim=c(-1,2), 
       border = FALSE)
  abline(v=a, col='#C79999', lwd=1)
  lines(density(priors$a), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$b, xlab="", prob=TRUE, 
       xlim=c(.5e3,2.5e3), 
       border = FALSE)
  abline(v=b, col='#C79999', lwd=1)
  lines(density(priors$b), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$c, xlab="", prob=TRUE, 
       xlim=c(0,1e3), 
       border = FALSE)
  abline(v=c, col='#C79999', lwd=1)
  lines(density(priors$c), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$sigma, xlab="", prob=TRUE, 
       xlim=c(-.1,.2), 
       border = FALSE)
  abline(v=sigma, col='#C79999', lwd=1)
  lines(density(priors$sigma), col='#6B8E8E', lwd=2, lty=3)
}

### question 1.2 ###

# 0. study design
# now, we have two populations, with different optimum GDDs
# + a kind of interaction, not sure about its biological signification though

gaussian2 <- function(x,g,a,b,db,c,d){
  return(0.7*exp(-(a*x-(b+g*db)+x*g*d)^2/(2*c^2)))
}

# 1. build the model, look at priors
model_file <- file.path(wd, 'model2.stan')
{
  priors <- list()
  priors$a <- rnorm(10000, 1, 0.5)
  priors$b <- rnorm(10000, 2000, 250)
  priors$db <- rnorm(10000, 0, 200)
  priors$c <- rnorm(10000, 500, 130)
  priors$d <- rnorm(10000, 0, 0.3)
  par(mfrow=c(1, 4))
  hist(priors$a, xlab="", prob=TRUE, 
       xlim=c(-1,2.5), 
       border = FALSE)
  hist(priors$b, xlab="", prob=TRUE, 
       xlim=c(1e3,3e3), 
       border = FALSE)
  hist(priors$db, xlab="", prob=TRUE, 
       xlim=c(-1e3,1e3), 
       border = FALSE)
  hist(priors$d, xlab="", prob=TRUE, 
       xlim=c(-1,1), 
       border = FALSE)
}

# 2. simulate data, check the model
a <- 0.8
b <- 1700 # optimal gdd
db <- 700 # higher optimal gdd for one group (g=1): 1300+700
c <- 600
d <- 0.4# 1/2*a
n = 100
simulate_data <- function(n, a, b, db, c, d, sigma = 0.05){
  x <- runif(n, 0, 5000) # simulate some annual GDD
  g <- sample(0:1, n, replace = TRUE)
  yhat <- 0.1 + gaussian2(x,g,a,b,db,c,d)
  y <- yhat + rnorm(n, mean = 0, sd = sigma)
  return(list(x=x,y=y,g=g))
}
sdata <- simulate_data(n, a, b, db, c, d)
par(mfrow=c(1, 1))
plot(sdata$y ~ sdata$x, col = c("#C79999","#6B8E8E")[sdata$g+1])

## what is the impact of d? it delays/advances the GDD optimum of pop. 1
# for(di in seq(-0.5,0.5,0.1)){
#   datad <- simulate_data(n, a, b, db, c, d = di)
#   par(mfrow=c(1, 1))
#   plot(datad$y ~ datad$x, col = c("#C79999","#6B8E8E")[datad$g+1])
# }


fit  <- stan(model_file, data=list(N=n, y= sdata$y, g=sdata$g, gdd = sdata$x), iter=2000, chains=4, cores = 4)
posteriors <- extract(fit, pars = c("a", "b", "db", "c", "d"))
{
  par(mfrow=c(1, 5))
  hist(posteriors$a, xlab="", prob=TRUE, 
       xlim=c(0,1.5), 
       border = FALSE)
  abline(v=a, col='#C79999', lwd=1)
  lines(density(priors$a), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$b, xlab="", prob=TRUE, 
       xlim=c(.5e3,2.5e3), 
       border = FALSE)
  abline(v=b, col='#C79999', lwd=1)
  lines(density(priors$b), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$db, xlab="", prob=TRUE, 
       xlim=c(-1e3,1e3), 
       border = FALSE)
  abline(v=db, col='#C79999', lwd=1)
  lines(density(priors$db), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$c, xlab="", prob=TRUE, 
       xlim=c(0,1e3), 
       border = FALSE)
  abline(v=c, col='#C79999', lwd=1)
  lines(density(priors$c), col='#6B8E8E', lwd=2, lty=3)
  hist(posteriors$d, xlab="", prob=TRUE, 
       xlim=c(-1,1), 
       border = FALSE)
  abline(v=d, col='#C79999', lwd=1)
  lines(density(priors$d), col='#6B8E8E', lwd=2, lty=3)
}

# question 2

# question 2.1, sampling from the data
sampling_size_effect <- 
  lapply(2:100, function(n){
    s <- sample(1:100,n,replace=FALSE)
    par(mfrow=c(1, 1))
    plot(sdata$y[s] ~ sdata$x[s], col = c("#C79999","#6B8E8E")[sdata$g[s]+1])
    fit  <- stan(model_file, data=list(N=n, y= sdata$y[s], g=sdata$g[s], gdd = sdata$x[s]), iter=2000, chains=4, cores = 4)
    posteriors <- extract(fit, pars = c("a", "b", "db", "c", "d"))
    data.frame(
      n = rep(n, 5),
      q5 = c(quantile(posteriors$a, 0.05), quantile(posteriors$b, 0.05), quantile(posteriors$db, 0.05), quantile(posteriors$c, 0.05), quantile(posteriors$d, 0.05)),
      q95 = c(quantile(posteriors$a, 0.95), quantile(posteriors$b, 0.95), quantile(posteriors$db, 0.95), quantile(posteriors$c, 0.95), quantile(posteriors$d, 0.95)),
      med = c(median(posteriors$a), median(posteriors$b),median(posteriors$db), median(posteriors$c), median(posteriors$d)),
      param = c('Parameter a', 'Parameter b', 'Parameter db', 'Parameter c', 'Parameter d')
    )
  })
sampling_size_effect <- as.data.frame(do.call(rbind, sampling_size_effect))

plot_question2.1 <- ggplot(data=sampling_size_effect) +
  facet_wrap(~ param, scales = 'free_x') +
  geom_pointrange(aes(y = n, xmin = q5, xmax = q95, x = med), 
                  col = 'grey70', size = 0.1) +
  geom_vline(data = data.frame(param = c('Parameter a', 'Parameter b', 'Parameter db', 'Parameter c', 'Parameter d'),
                               x = c(a,b, db,c,d)),
             aes(xintercept = x), col = '#C79999') +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), axis.line = element_line())
ggsave(plot_question2.1, filename = file.path(wd, 'plot_question2.1.pdf'), height = 6, width = 10)

# question 2.2, repeated sampling from the data
sampling_size_effect2 <- 
  lapply(seq(5,85,10), function(n){
    posteriors <- lapply(1:5, function(i){ # only 5 samples because it's tooooo long
      s <- sample(1:100,n,replace=FALSE)
      par(mfrow=c(1, 1))
      plot(sdata$y[s] ~ sdata$x[s], col = c("#C79999","#6B8E8E")[sdata$g[s]+1])
      fit  <- stan(model_file, data=list(N=n, y= sdata$y[s], g=sdata$g[s], gdd = sdata$x[s]), iter=2000, chains=4, , cores = 4)
      posteriors <- extract(fit, pars = c("a", "b", "db", "c", "d"))
      data.frame(posteriors)
    })
    posteriors <- as.data.frame(do.call(rbind, posteriors))
    
    data.frame(
      n = rep(n, 5),
      q5 = c(quantile(posteriors$a, 0.05), quantile(posteriors$b, 0.05), quantile(posteriors$db, 0.05), quantile(posteriors$c, 0.05), quantile(posteriors$d, 0.05)),
      q95 = c(quantile(posteriors$a, 0.95), quantile(posteriors$b, 0.95), quantile(posteriors$db, 0.95), quantile(posteriors$c, 0.95), quantile(posteriors$d, 0.95)),
      med = c(median(posteriors$a), median(posteriors$b),median(posteriors$db), median(posteriors$c), median(posteriors$d)),
      param = c('Parameter a', 'Parameter b', 'Parameter db', 'Parameter c', 'Parameter d')
    )
  })
sampling_size_effect2 <- as.data.frame(do.call(rbind, sampling_size_effect2))

plot_question2.2 <- ggplot(data=sampling_size_effect2) +
  facet_wrap(~ param, scales = 'free_x') +
  geom_pointrange(aes(y = n, xmin = q5, xmax = q95, x = med), 
                  col = 'grey70', size = 0.1) +
  geom_vline(data = data.frame(param = c('Parameter a', 'Parameter b', 'Parameter db', 'Parameter c', 'Parameter d'),
                               x = c(a,b, db,c,d)),
             aes(xintercept = x), col = '#C79999') +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), axis.line = element_line())
ggsave(plot_question2.2, filename = file.path(wd, 'plot_question2.2.pdf'), height = 6, width = 10)

# question 3
n <- 100
sdata <- simulate_data(n, a, b, db, c, d, sigma = 0.05*1.5) # increase error term by 50%
sdata$y[sdata$y<0] <- 0 # fitness cannot be negative?

sampling_size_effect_error <- 
  lapply(seq(5,85,10), function(n){ 
    posteriors <- lapply(1:5, function(i){ # only 5 samples because it's tooooo long
      s <- sample(1:100,n,replace=FALSE)
      par(mfrow=c(1, 1))
      plot(sdata$y[s] ~ sdata$x[s], col = c("#C79999","#6B8E8E")[sdata$g[s]+1])
      fit  <- stan(model_file, data=list(N=n, y= sdata$y[s], g=sdata$g[s], gdd = sdata$x[s]), iter=2000, chains=4, cores = 4)
      posteriors <- extract(fit, pars = c("a", "b", "db", "c", "d"))
      data.frame(posteriors)
    })
    posteriors <- as.data.frame(do.call(rbind, posteriors))
    
    data.frame(
      n = rep(n, 5),
      q5 = c(quantile(posteriors$a, 0.05), quantile(posteriors$b, 0.05), quantile(posteriors$db, 0.05), quantile(posteriors$c, 0.05), quantile(posteriors$d, 0.05)),
      q95 = c(quantile(posteriors$a, 0.95), quantile(posteriors$b, 0.95), quantile(posteriors$db, 0.95), quantile(posteriors$c, 0.95), quantile(posteriors$d, 0.95)),
      med = c(median(posteriors$a), median(posteriors$b),median(posteriors$db), median(posteriors$c), median(posteriors$d)),
      param = c('Parameter a', 'Parameter b', 'Parameter db', 'Parameter c', 'Parameter d')
    )
  })
sampling_size_effect_error <- as.data.frame(do.call(rbind, sampling_size_effect_error))

plot_question3 <- ggplot(data=sampling_size_effect_error) +
  facet_wrap(~ param, scales = 'free_x') +
  geom_pointrange(aes(y = n, xmin = q5, xmax = q95, x = med), 
                  col = 'grey70', size = 0.1) +
  geom_vline(data = data.frame(param = c('Parameter a', 'Parameter b', 'Parameter db', 'Parameter c', 'Parameter d'),
                               x = c(a,b, db,c,d)),
             aes(xintercept = x), col = '#C79999') +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), axis.line = element_line())
ggsave(plot_question3, filename = file.path(wd, 'plot_question3.pdf'), height = 6, width = 10)
