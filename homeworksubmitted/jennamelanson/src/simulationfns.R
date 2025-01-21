######################################
#### Functions for data simulation
######################################

simulateLinRegData <- function(sample_size, 
                               intercept, 
                               coefs, 
                               predictors, 
                               sigma,
                               interaction_coefs = NULL){
  ##IMPORTANT: COEFS and PREDICTORS must be in the same order, and the same length
  
  if (length(coefs) != length(predictors)) {
    stop("Error: the number of coefficients does not match the number of predictors.")
  }
  
  #initialize ypred and output df
  ypred = rep(intercept, sample_size)
  df = data.frame(rep(NA,sample_size))
  
  #run main effects
  for (i in 1:length(predictors)) {
    ypred = ypred + coefs[i]*predictors[[i]]
    
    cname = paste("pred", i, sep = "")
    df$V1 = predictors[[i]]
    names(df)[names(df)=="V1"] = cname
  }
  
  #run all two way interactions
  # the correct order to input interaction coefficients is:
  #for four predictors: 1x2, 1x3, 1x4, 2x3, 2x4, 3x4
  if(!is.null(interaction_coefs)){
    init = 1
    for (j in 1:length(predictors)){
      k = j+1
        while (k <= length(predictors)) {
          ypred = ypred + interaction_coefs[init]*predictors[[j]]*predictors[[k]]
          k <- k + 1
          init = init + 1
        }
      }
    }
  
  #calculate yobs
  yobs = rnorm(sample_size, ypred, sigma)
  
  #create dataframe with predictors and response as output
  df$yobs = yobs
  df = df[,2:length(colnames(df))]
  return(df)
}



simulateChallenge1 <- function(sample_size, 
                               intercept,
                               sigma,
                               b1,
                               btreatment1,
                               btreatment2,
                               interaction_treatment1,
                               interaction_treatment2){
 #This function is NOT very generalizable, but I am simply too lazy to make it so
  
  
  #next, generate predictor 1--uniform distribution
  pred1 = runif(sample_size, 0, 1)
  
  #next, generate data with 3 levels
  treatment = c("control", "treatment1", "treatment2")
  pred2 = as.data.frame(treatment) %>% sample_n(sample_size, replace = TRUE)
  if(length(unique(pred2$treatment)) < 3){
    pred2$treatment[1] = "control"
    pred2$treatment[2] = "treatment1"
    pred2$treatment[3] = "treatment2"
  }
  
  #join into df
  predictor_df = data.frame(pred1, pred2)
  print(predictor_df$treatment)
  
  #use model_matrix to convert possible_treatments into two columns of 0/1
  dummy <- model.matrix(~factor(treatment), data = predictor_df)
  predictor_df$treatment1 = dummy[,2]
  predictor_df$treatment2 = dummy[,3]
  
  #get ypred
  ypred = intercept + 
    b1*predictor_df$pred1 + 
    btreatment1*predictor_df$treatment1 + 
    btreatment2*predictor_df$treatment2 +
    interaction_treatment1*pred1*predictor_df$treatment1 +
    interaction_treatment2*pred1*predictor_df$treatment2
  #the third level (control) has no effect, so doesn't even need to be coded in here
  
  #simulate yobs
  yobs = rnorm(sample_size, ypred, sigma)
  data = data.frame(predictor_df, yobs)
  
  return(data)
}




simulateChallenge2 <- function(sample_size, 
                               intercept,
                               sigma,
                               beta){
  #This function is NOT very generalizable, but I am simply too lazy to make it so
  
  
  #generate predictor 1--uniform distribution
  pred1 = runif(sample_size, 0, 1)
  
  #next, generate two binary variables -- bernoulli distribution
  pred2 = rbinom(sample_size, 1, prob = 0.5)
  pred3 = rbinom(sample_size, 1, prob = 0.5)
  
  #join into df
  predictor_df = data.frame(pred1, pred2, pred3)
  
  #get ypred
  ypred = intercept + 
    beta[1]*predictor_df$pred1 + 
    beta[2]*predictor_df$pred2 + 
    beta[3]*predictor_df$pred3 +
    beta[4]*predictor_df$pred1*predictor_df$pred2 + #two way inter (1vs2)
    beta[5]*predictor_df$pred1*predictor_df$pred3 + #two way inter (1vs3)
    beta[6]*predictor_df$pred2*predictor_df$pred3 + #two way inter (2vs3)
    beta[7]*predictor_df$pred1*predictor_df$pred2*predictor_df$pred3 #three way inter (1vs2vs3)
  
  #simulate yobs
  yobs = rnorm(sample_size, ypred, sigma)
  data = data.frame(predictor_df, yobs)
  
  return(data)
}


######################################
#### Functions for monte carlo
######################################

monteCarloChallenge1 <- function(number_repeats,
                       sample_size, 
                       intercept,
                       sigma,
                       b1 = NULL,
                       btreatment1 = NULL,
                       btreatment2 = NULL,
                       interaction_treatment1 = NULL,
                       interaction_treatment2 = NULL
                       ){
  #define sample size range
  n = 3:sample_size
  
  #initialize dataframe to store parameter estimates
  df_size = number_repeats*(sample_size - 2)
  param_df = data.frame(rep(NA,df_size),rep(NA,df_size),rep(NA,df_size),rep(NA,df_size), rep(NA,df_size), rep(NA,df_size),rep(NA,df_size),rep(NA,df_size))
  colnames(param_df) = c("sample_size", "sigma", "intercept", "b1", "btreatment1", "btreatment2", "interaction_treatment1", "interaction_treatment2")
  
  #initialize count (for recording results)
  count = 1
  
  #simulate + fit models!
  for (i in n){
    for (j in 1:number_repeats){
      
      #simulate data
      simData = simulateChallenge1(sample_size = i,
                                   intercept = intercept,
                                   sigma = sigma,
                                   b1 = b1,
                                   btreatment1 = btreatment1,
                                   btreatment2 = btreatment2,
                                   interaction_treatment1 = interaction_treatment1,
                                   interaction_treatment2 = interaction_treatment2)
      
      #check that predictors aren't constant
      #make sure each treatment level occurs at least once (coerce if not)
      if(length(unique(simData$treatment)) < 3){
        simData$treatment[1] = "control"
        simData$treatment1[1] = 0
        simData$treatment2[1] = 0
        simData$treatment[2] = "treatment1"
        simData$treatment1[2] = 1
        simData$treatment2[2] = 0
        simData$treatment[3] = "treatment2"
        simData$treatment1[3] = 0
        simData$treatment2[3] = 1
      }
      
      #fit model
      fit = rstanarm::stan_glm(yobs ~ pred1*treatment1 + pred1*treatment2, data = simData)
      
      #extract parameter estimates
      param_df[count,1] = i
      param_df[count, 2] = sigma(fit)
      param_df[count, 3] = coef(fit)[1]
      param_df[count, 4] = coef(fit)[2]
      param_df[count, 5] = coef(fit)[3]
      param_df[count, 6] = coef(fit)[4]
      param_df[count, 7] = coef(fit)[5]
      param_df[count, 8] = coef(fit)[6]
      count = count + 1
    }
  }
  
  return(param_df)
}




monteCarloChallenge2 <- function(number_repeats,
                                 sample_size, 
                                 intercept,
                                 sigma,
                                 beta,
                                 confidence_interval
){
  #define sample size range
  #start at four because we want a fully factorial design
  n = 4:sample_size
  
  #initialize dataframe to store parameter estimates
  df_rows = number_repeats*(sample_size - 3)
  number_params = 2 + length(beta)
  df_cols = 1 + number_params*3
  param_df = data.frame((matrix(ncol = df_cols, nrow = df_rows)))
  colnames(param_df) = c("sample_size", 
                         "sigma", "sigmalower", "sigmaupper",
                         "intercept", "interceptlower", "interceptupper",
                         "b1", "b1lower", "b1upper",
                         "b2", "b2lower", "b2upper",
                         "b3", "b3lower", "b3upper",
                         "b4", "b4lower", "b4upper",
                         "b5", "b5lower", "b5upper",
                         "b6", "b6lower", "b6upper",
                         "b7", "b7lower", "b7upper")
  
  #initialize count (for recording results)
  count = 1
  
  #simulate + fit models!
  for (i in n){
    for (j in 1:number_repeats){
      
      #simulate data
      simData = simulateChallenge2(sample_size = i,
                                   intercept = intercept,
                                   sigma = sigma,
                                   beta = beta)
      
      #check that predictors aren't constant for small sample size
      #coerce to fully factorial
      if(length(unique(paste(simData$pred2, simData$pred3)) < 4)){
        simData$pred2[1] = 1
        simData$pred3[1] = 0
        simData$pred2[2] = 0
        simData$pred3[2] = 1
        simData$pred2[3] = 0
        simData$pred3[3] = 0
        simData$pred2[4] = 1
        simData$pred3[4] = 1
      }
      
      #fit model
      fit = rstanarm::stan_glm(yobs ~ pred1*pred2*pred3, data = simData)
      intervals = posterior_interval(fit, prob = confidence_interval)
      
      #extract parameter estimates
      #surely there must be a more efficient way to do this lol
      param_df[count,1] = i
      param_df[count, 2] = sigma(fit)
      param_df[count,3:4] = intervals[rownames(intervals) == "sigma"]
      param_df[count, 5] = coef(fit)[1]
      param_df[count,6:7] = intervals[rownames(intervals) == "(Intercept)"]
      param_df[count, 8] = coef(fit)[2]
      param_df[count,9:10] = intervals[rownames(intervals) == "pred1"]
      param_df[count, 11] = coef(fit)[3]
      param_df[count,12:13] = intervals[rownames(intervals) == "pred2"]
      param_df[count, 14] = coef(fit)[4]
      param_df[count,15:16] = intervals[rownames(intervals) == "pred3"]
      param_df[count, 17] = coef(fit)[5]
      param_df[count,18:19] = intervals[rownames(intervals) == "pred1:pred2"]
      param_df[count, 20] = coef(fit)[6]
      param_df[count,21:22] = intervals[rownames(intervals) == "pred1:pred3"]
      param_df[count, 23] = coef(fit)[7]
      param_df[count,24:25] = intervals[rownames(intervals) == "pred2:pred3"]
      param_df[count, 26] = coef(fit)[8] 
      param_df[count,27:28] = intervals[rownames(intervals) == "pred1:pred2:pred3"]
      count = count + 1
      
    }
  }
  
  return(param_df)
}