set.seed(1982)
library(R2jags)
library(polspline)

# This script applies the CC, EWA, and LMM models
# to data from the paper Hunger affects social decisions...
# Pre-print available https://psyarxiv.com/67abq/
# Data available 
# https://figshare.com/articles/Data_and_code_for_Fraser_and_Nettle_Hunger_affects_social_decisions_in_a_Public_Goods_Game_but_not_an_Ultimatum_Game_/11807346/1

# analysis is not conducted on hunger condition. Data are used for model comparison

nsubs <- 4
groupSize <- nsubs
ntrials <- 10
pi <- 1.4
ntokens <- 20
vals <- seq(0,20,1) #possible values to contribute - from 0 to 20 tokens

setwd("SET TO DATA DIRECTORY")

rawDat <- read.csv("hungerAffectsSocialDecisions.csv") # Public goods game

group_names <- unique(rawDat$UniqueGroup)
ngroups <- length(group_names)

#---------------- create data matrices - subject x trial x group x condition-------------------------------------------------

# data for no punishment condition #
c_no_punish <- array(0,c(nsubs,ntrials,ngroups))
Ga_no_punish <- array(0,c(ntrials,ngroups))
Gc_no_punish <- array(0,c(nsubs,ntrials,ngroups))

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="N"][1:10],
                            rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="N"][11:20],
                            rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="N"][21:30],
                            rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="N"][31:40])
  
  Ga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  for (s in 1:nsubs) {
    Gc_no_punish[,,g] <- colSums(c_no_punish[-s,,g])
  }
}

# desired data matrices
c <- c_no_punish
Ga <- Ga_no_punish
Gc <- Gc_no_punish
c_choice_index <- c

#------------ Model comparison --------------------------------------------------------------------------

# apply hierarchical versions of CC, EWA, LMM, and motivational EWA models 

setwd("SET TO MODEL FILE DIRECTORY")

data <- list("groupSize", "ngroups", "ntrials", "ntokens", "pi", "vals","c","Gc","c_choice_index","Ga") #data inputted into jags
params <- c("mu_c") #parameters we'll track in jags
samplesCC <- jags(data, inits=NULL, params,
                model.file ="CC_group.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

params <- c("mu_c") #parameters we'll track in jags
samplesEWA <- jags(data, inits=NULL, params,
                  model.file ="EWA_group.txt",
                  n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

params <- c("nu_c") #parameters we'll track in jags
samplesLMM <- jags(data, inits=NULL, params,
                   model.file ="LMM_group.txt",
                   n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

# Calculate rewards seperately, to apply to utility function for EWA model
R_self <- ntokens - c + ((pi/groupSize)*(Gc+c)) # calculate own outcome
R_other <- array(0,c(nsubs,ntrials,ngroups)) #calculate outcome for others, assuming only knowledge of average group contrib.  
R_other[1,,] <- ntokens - Ga + ((pi/groupSize)*(Gc[1,,]+Ga)) #subject 1  
R_other[2,,] <- ntokens - Ga + ((pi/groupSize)*(Gc[2,,]+Ga)) #subject 2  
R_other[3,,] <- ntokens - Ga + ((pi/groupSize)*(Gc[3,,]+Ga)) #subject 3  
R_other[4,,] <- ntokens - Ga + ((pi/groupSize)*(Gc[4,,]+Ga)) #subject 4  

data <- list("groupSize", "ngroups", "ntrials", "ntokens", "pi", "vals","c","Gc","c_choice_index","Ga","R_self","R_other") #data inputted into jags
params <- c("mu_c") #parameters we'll track in jags
samplesEWAmotives <- jags(data, inits=NULL, params,
                          model.file ="EWA_motives_group.txt",
                          n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)
