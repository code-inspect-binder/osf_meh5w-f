LMM_fun <- function(nagents,ntrials,vals,ntokens,pi,parameters) { 
  
  #Z <- c(1,1,0) #indicates whether agent is conditional co-operator or not
  
  #---------------- parameters ---------------------------
  #CC model
  Gb1 <- parameters$Gb1
  omega1 <- parameters$omega1 #initial weighting of beliefs about others contributions in choice of own contribution, relative to prefs
  lambda <- parameters$lambda #decay rate in weighting of beliefs about others - prefs/predictions dominate over time
  gamma <- parameters$gamma #parameter weighting of beliefs about what others will contribute, relative to observed contribution
  p0 <- parameters$p0 #intercept of linear model relating preferred contributions to possible contribution values
  pbeta <- parameters$pbeta #slope of linear model relating preferred contributions to possible contribution values

  #EWA model
  delta <- parameters$delta # weighting of forgone versus received payoffs
  rho <- parameters$rho #higher number means older trials have great influence and strategy fixation is slower
  phi <- parameters$phi #memory of old attractions - volatility assumption about environ. priors (0,1)
  theta <- parameters$theta   #consistency of choice with attractions - inverse heat - explore exploit
  
  Z <- parameters$Z
  
  #-----------simulation arrays - to be filled------------ 
  #CC model
  Ga <- array(0,c(ntrials)) #observed others' contribution (mean of all others)
  Gb <- array(0,c(nagents,ntrials)) #beliefs about others' contribution (mean of all others)
  p <- array(0,c(nagents,ntrials)) #preffered contribution on each trial (independent of beliefs)
  omega <- array(0,c(nagents,ntrials)) #weighting of beliefs about others' contribution relative to prefs
  c_CC <- array(0,c(nagents,ntrials)) #actual contributions

  #EWA model
  N <- array(0,c(nagents,ntrials))
  A <- array(0,c(nagents,ntrials,ntokens+1))
  expA <- array(0,c(nagents,ntrials,ntokens+1))
  P <- array(0,c(nagents,ntrials,ntokens+1))
  c_EWA_index <- array(0,c(nagents,ntrials))
  c_EWA <- array(0,c(nagents,ntrials))
  
  c <- array(0,c(nagents,ntrials))
  
  #-----------initialisations - first trials------------ 
  # CC model
  # agents preferences - assumed to be a linear function of possible values - linear function has p0 and pbeta as params
  pvals <- array(0,c(nagents,length(vals)))
  for (n in 1:nagents) {
    pvals[n,] <- p0[n] + (pbeta[n]*vals) #vector of preferred contributions for each possible value - assume linear relationship
    # this deviates from paper, which also has "triangular" realationships for some (i.e. increase co-operation up to a value then <)
  }
  # set omega as starting weighting for beliefs relative to preferences as parameter
  omega[,1] <- omega1
  # set starting beliefs about what others will contribute
  Gb[,1] <- Gb1
  # set average first trial contribution to average belief, assume full co-operation at outset. 
  # Reasonable simplification.
  c_CC[,1] <- Gb1
  Ga[1] <- mean(Gb1)
  
  # EWA mode
  N[,1] <- 1
  A[,1,] <- 0#rnorm(10,10,10)
  c_EWA[,1] <- 20
  
  c[,1] <- 20
  
  #-----------simulation----------------------------------- 

  for (t in 2:ntrials) {
    
    ################################################################################
    #################### CC model ##################################################
    ################################################################################
    
    for (n in 1:nagents) {
      
      # update beliefs about what others contribute as a weighted average of prior beliefs and observed contributions
      #from page 549 "subjects belief in a given period is a weighted average of what he or she believed 
      #about others in the previous period and his or her observations of others' contributions" 
      Gb[n,t] <- (gamma[n]*(Gb[n,t-1]))+((1-gamma[n])*(Ga[t-1]))
      
      # determine what people "predict" or "prefer" to contribute, given the group contribution, outside
      # of, or independent of, their willingness to co-operate. Index preferences relative to believed contribution of others
      p[n,t] <- pvals[n,round(Gb[n,t])]
      
      # update relative weighting of beliefs about others contributions, using decay function
      #page 550 "in later periods, predicted (i.e. preferred) contibution becomes more important than belief"
      omega[n,t] <- omega[n,t-1]*(1-lambda[n]) #1 - lamba - because lambda is a decay rate
      # departure from paper - possible model innovation
      
      #page 550 subjects contribute a weighted average of predicted (i.e. preferred) contribution and belief
      c_CC[n,t] <- ceiling((omega[n,t])*Gb[n,t] + ((1-omega[n,t])*p[n,t])) 
      
    }
    
    
    ################################################################################
    #################### EWA model ##################################################
    ################################################################################

    for (n in 1:nagents) {
      
      N[n,t] <- (rho[n]*N[n,t-1]) + 1 #eqn 2.1
      
      for (j in 1:(ntokens+1)) { # + 1 because there are n + 1 possible contribs - integer vals & 0 - need to treat A as index 
        
        A[n,t,j] <- ( 
          
          (phi[n]*N[n,t-1]*A[n,t-1,j]) + #prior attractions
            (delta[n] + ((1-delta[n])*(c[n,t-1]==j))) * #indicates whether jth token was chosen
            ((((j+sum(c[-n,t-1]))*pi)/nagents)-j) # calculate payoff for each possible contrib.
          
        )/
          N[n,t] #experience weighting
        
        expA[n,t,j] <- exp(theta[n]*A[n,t,j])
        
      }
      
      for (j in 1:(ntokens+1)) {    
        P[n,t,j] <- expA[n,t,j]/sum(expA[n,t,])
      }
      
      c_EWA_index[n,t] <- rcat(1,P[n,t,])
      
      c_EWA[n,t] <- vals[c_EWA_index[n,t]]
      
    }
    
    #---------------------assign contribution depending on model ----------------------------
    for (n in 1:nagents) {
      c[n,t] <- ifelse(Z[n]==1,c_CC[n,t],c_EWA[n,t])
    }
    
    #---------------------averages and beliefs for models -----------------------------------
    #CC model
    # recode average contribution as observed belief about how much each agent contributed Ga
    Ga[t] <- sum(c[,t])/nagents
  
    # # calculate *others* contributions from contribution matrix c, to enter as data
    # # each row vector for agent n represents the average of what the *others* contributed on the trial
    # # used in atraction equation in jags instead of sum(c[-n],t-1), because no negative indexing in jags
    # Gc <- array(0,c(nagents,ntrials))
    # for (n in 1:nagents) {
    #   Gc[n,] <- colSums(c[-n,])
    # }
    
  }
  
  result <- list(c=c,c_CC=c_CC,c_EWA,c_EWA)
  
  return(result)
  
  
}
  
  