set.seed(1982)
library(R2jags)
library(polspline)

groupSize <- 4
ntrials <- 10
pi <- 1.4
ntokens <- 20
vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

setwd("SET DATA DIRECTORY")

rawDat <- read.csv("AntisocialPubishmentAcrossSocieties.csv") # Public goods game

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

groupSize <- 4
ntrials <- 10
pi <- 1.4
ntokens <- 20
vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

setwd("SET MODEL FILE DIRECTORY")

rawDat <- read.csv("AntisocialPubishmentAcrossSocieties.csv") # Public goods game

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups))
Ga_no_punish <- array(0,c(ntrials,ngroups))
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups))
missing <- array(0,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Ga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[,,g] <- colSums(c_no_punish[-s,,g])
  }
}

# data for punishment condition #
c_punish <- array(0,c(groupSize,ntrials,ngroups))
Ga_punish <- array(0,c(ntrials,ngroups))
Gc_punish <- array(0,c(groupSize,ntrials,ngroups))

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Ga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    Gc_punish[,,g] <- colSums(c_punish[-s,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Ga <- array(0,c(ntrials,ngroups,2))
Ga[,,1] <- Ga_no_punish
Ga[,,2] <- Ga_punish

Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

c_choice_index <- c

civic <- array(0,c(ngroups))
for (g in 1:ngroups) {
  civic[g] <- mean(redDat$civic[redDat$groupid==group_names[g]&redDat$p=="P-experiment"])
}

civic <- array(0,c(ngroups))
for (g in 1:ngroups) {
  civic[g] <- mean(redDat$civic[redDat$groupid==group_names[g]&redDat$p=="P-experiment"])
}

# remove data from coutries which do not have civic scores
c <- c[,,!is.na(civic),]
Ga <- Ga[,!is.na(civic),]
Gc <- Gc[,,!is.na(civic),]
civic <- civic[!is.na(civic)]

#standardise civic norms score
civic <- (civic-mean(civic))/sd(civic)

#redefine number of groups after removing those without civic scores
ngroups <- length(civic)

data <- list("groupSize", "ngroups", "ntrials", "ntokens", "pi", "vals","c","Gc","c_choice_index","Ga","civic") #data inputted into jags
params <- c("alpha","mu","beta0","betaC") #parameters we'll track in jags

# - run jags code
samples <- jags(data, inits=NULL, params,
                model.file ="civic_norms_GLM.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)


prior <- dnorm(0,1)
fit.posterior <- logspline(samples$BUGSoutput$sims.list$betaC)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior/posterior

qlogspline(0.025,fit.posterior) #2.5% CI
qlogspline(0.975,fit.posterior) #97.5% CI


prior <- dnorm(0,1)
fit.posterior <- logspline(samples$BUGSoutput$sims.list$beta0)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior/posterior

qlogspline(0.025,fit.posterior) #2.5% CI
qlogspline(0.975,fit.posterior) #97.5% CI

par(mfrow=c(1,3),cex=1.5)
plot(density(samples$BUGSoutput$sims.list$betaC),xlab = expression(paste(beta["C"])),ylab=" ",
     #cex.lab=2, cex.axis = 2,
     lwd = 3,
     main=expression(paste("Posterior for ", beta["C"])),frame=FALSE) 
plot(civic,colMeans(colMeans(pnorm(samples$BUGSoutput$sims.list$mu))),ylim=c(0,1),
     #cex.lab=2, cex.axis = 2, cex.main = 2,
     xlab="Adherence to Civic Norms",
     ylab="",
     main=expression(paste("Prob. Coop.  -  Standard  -  ", mu)),frame=FALSE)
plot(civic,colMeans(colMeans(pnorm(samples$BUGSoutput$sims.list$mu+samples$BUGSoutput$sims.list$alpha))),ylim=c(0,1),
     #cex.lab=2, cex.axis = 2,cex.main = 2,
     ylab="",
     xlab="Adherence to Civic Norms",
     main=expression(paste("Prob. Coop.  -  Punish  -  ", mu, " + ", alpha)),frame=FALSE)



