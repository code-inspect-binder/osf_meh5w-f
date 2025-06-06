set.seed(1982)
library(extraDistr)
library(R2jags)

niterations <- 100
ntrials <- 10
nagents <- 10

# simulation information 
ntokens <- 20
pi <- 1.5
vals <- seq(0,20,1) #possible values to contribute - from 0 to 20 tokens

true <- c()
true$omega1 <- array(0,c(niterations,nagents))
true$lambda <- array(0,c(niterations,nagents))
true$gamma <- array(0,c(niterations,nagents))
true$pbeta <- array(0,c(niterations,nagents))

true$delta <- array(0,c(niterations,nagents))
true$rho <- array(0,c(niterations,nagents))
true$phi <- array(0,c(niterations,nagents))
true$theta <- array(0,c(niterations,nagents))

true$Z <- array(0,c(niterations,nagents))


infer <- c()
infer$omega1 <- array(0,c(niterations,nagents))
infer$lambda <- array(0,c(niterations,nagents))
infer$gamma <- array(0,c(niterations,nagents))
infer$pbeta <- array(0,c(niterations,nagents))

infer$delta <- array(0,c(niterations,nagents))
infer$rho <- array(0,c(niterations,nagents))
infer$phi <- array(0,c(niterations,nagents))
infer$theta <- array(0,c(niterations,nagents))

psy <- array(0,c(niterations,nagents))
Z <- array(0,c(niterations,nagents))

for (i in 1:niterations) { 

  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------- Set model parameters -------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  #-------------- CC model parameters ------------------------
  parameters <- c()
  #initial beliefs about others average contribution. Set to max as simplification. Could have been free/estimated
  parameters$Gb1 <- rep(20,nagents) 
  #initial weighting of beliefs about others contributionsrelative to own prefs. Higher number means > social influence
  parameters$omega1 <- runif(nagents,.9,1) 
  #decay rate in weighting of beliefs about others - prefs dominate over time following decay function
  parameters$lambda <- rbeta(nagents,1,1)#runif(nagents,.7,1)
  #parameter weighting of beliefs about what others will contribute, relative to observed contribution (learning rate)
  parameters$gamma <- rbeta(nagents,1,1)#runif(nagents,.5,1) 
  #intercept of linear model relating preferred contributions to possible contribution values
  parameters$p0 <- rep(0,nagents) 
  #slope of linear model relating preferred contributions to possible contribution values
  #this is capped at .7, in parameter recovery. values higher than this leave no room for beliefs because prefs are so high
  #should be free parameter in inference though
  parameters$pbeta <- rbeta(nagents,1,1)#runif(nagents,.5,1)#runif(nagents,.5,.7) 
  
  #-------------- EWA model parameters ------------------------
  #set free parameters - we can do inference on these
  parameters$delta <- rbeta(nagents,1,1)#runif(nagents,.1,1) 
  parameters$rho <- rbeta(nagents,1,1)#runif(nagents,.1,1)
  parameters$phi <- rbeta(nagents,1,1)#runif(nagents,.1,1)
  parameters$theta <- rgamma(nagents,1,1)#runif(3,.1,2) 
  #-----------------------------------------------------------
  
  # mixture parameter
  parameters$Z <- rbinom(10,1,.5)
  
  true$omega1[i,] <- parameters$omega1
  true$lambda[i,] <- parameters$lambda
  true$gamma[i,] <- parameters$gamma
  true$pbeta[i,] <- parameters$pbeta
  
  true$delta[i,] <- parameters$delta
  true$rho[i,] <- parameters$rho
  true$phi[i,] <- parameters$phi
  true$theta[i,] <- parameters$theta

  true$Z[i,] <- parameters$Z
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------- Run simulation -------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  # load simulation function into workspace
  setwd("SET MODEL FILE DIRECTORY")
  source("LMM_fun.R")
  sims <- LMM_fun(nagents,ntrials,vals,ntokens,pi,parameters)
  
  c <- sims$c
  
  #------ data for CC model inference
  Ga <- colMeans(c) #average group member contribution
  
  #-------data for EWA model inference
  # calculate *others* contributions from contribution matrix c, to enter as data
  # each row vector for agent n represents the average of what the *others* contributed on the trial
  # used in atraction equation in jags instead of sum(c[-n],t-1), because no negative indexing in jags
  Gc <- array(0,c(nagents,ntrials))
  for (n in 1:nagents) {
    Gc[n,] <- colSums(c[-n,])
  }
  
  # re-representation of contribution not as data but as model input, for attraction equation
  c_choice_index <- c
  
  data <- list("nagents", "ntrials", "ntokens", "pi", "vals","c","Gc","c_choice_index","Ga") #data inputted into jags
  params <- c("omega1", "lambda", "gamma", "pbeta", "delta", "rho", "phi", "theta", "Z","psy") #parameters we'll track in jags
  
  samples <- jags(data, inits=NULL, params,
                  model.file ="LMM.txt",
                  n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)
  
  # save maximum a posteriori (MAP) values for parameters from fitted model
  for (n in 1:nagents) {

    #-------------- CC model--------------------------
    X <- samples$BUGSoutput$sims.list$omega1[,n]
    infer$omega1[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$lambda[,n]
    infer$lambda[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$gamma[,n]
    infer$gamma[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$pbeta[,n]
    infer$pbeta[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    #-------------- EWA mode---------------------------    
    X <- samples$BUGSoutput$sims.list$delta[,n]
    infer$delta[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$rho[,n]
    infer$rho[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$phi[,n]
    infer$phi[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$theta[,n]
    infer$theta[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]

    X <- samples$BUGSoutput$sims.list$psy[,n]
    psy[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))]
    
    X <- samples$BUGSoutput$sims.list$Z[,n]
    Z[i,n] <-density(X)$x[which(density(X)$y==max(density(X)$y))][1]
    
  }
  
  print(i)
  
}

accuracy <- sum(round(Z)==true$Z)
HR <- sum(round(Z)==1&true$Z==1)/sum(true$Z==1)
FA <- sum(round(Z)==1&true$Z==0)/sum(true$Z==0)