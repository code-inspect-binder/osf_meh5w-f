set.seed(1982)
library(R2jags)
library(polspline)

# This script applies a hierarchical t-test over latent mixture model parameters
# to data from the paper Prosocial preferences do no explain human cooperation in public goods games
# Paper available 
# https://www.pnas.org/content/110/1/216
# Companion paper
# https://royalsocietypublishing.org/doi/full/10.1098/rspb.2014.2678
# Data available 
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.cr829

# Design is odd - co-operating groups are not recurring but are reformed from other members
# in a session on each trial - expect generally lower conditional co-operation
# Analysis is not conducted on learning condition/expectatons/blackbox condition. 
# Data are used to analyse only effects of different MPCR values on standard public goods game

nsubs <- 4
groupSize <- nsubs
ntrials <- 20
pi <- c(1.6,6.4)
ntokens <- 40
vals <- seq(1,41,1) #possible values to contribute - from 0 to 40 tokens

setwd("SET DATA DIRECTORY")

rawDat <- read.csv("LearningInformation.csv") # Public goods game

rawDat$sumOthersContrib <- (rawDat$otherA+rawDat$otherB+rawDat$otherC) #calculate others contribution, will need

# in this experiment, participants are not put in stable groups. instead, groups re-sort during each session
# required to input data by subject and not by group, since no real groups exist
# but data is recorded for total group contribution, so all info required for models is available
# should expect different social dynamics 

# in this analysis, we compare conditional cooperation between two levels of MPCR in the experiment - the
# MCPR = .4 (multiplication factor pi = 1.6) and MCPR = 1.6 (multiplication factor pi = 6.4) conditions - in the standard 
# public goods game. Design is within subjects

# extract data for the pi = 1.6 multiplication factor trials
rawDat_16 <- rawDat[rawDat$EF==1.6 & 
                          rawDat$Game=="rPG" & 
                          rawDat$Treatment=="Partial",]

# extract data for the pi = 6.4 multiplication factor trials
rawDat_64 <- rawDat[rawDat$EF==6.4 & 
                      rawDat$Game=="rPG" & 
                      rawDat$Treatment=="Partial",]

rawDat_16$uniqueSubno <- rep(1:(length(rawDat_16$EF)/20),each = 20)
rawDat_64$uniqueSubno <- rep(1:(length(rawDat_64$EF)/20),each = 20)

n_16 <- max(rawDat_16$uniqueSubno)
n_64 <- max(rawDat_64$uniqueSubno)

# calculate data arrays for pi = 1.6
c_16 <- c()
Ga_16 <- c()
Gc_16 <- c()
for (i in 1:n_16) { 
  c_16 <- rbind(c_16,rawDat_16$Contribution[rawDat_16$uniqueSubno==i])
  Ga_16 <- rbind(Ga_16,rawDat_16$averageC[rawDat_16$uniqueSubno==i])
  Gc_16 <- rbind(Gc_16,rawDat_16$sumOthersContrib[rawDat_16$uniqueSubno==i])
}

# calculate data arrays for pi = 6.4
c_64 <- c()
Ga_64 <- c()
Gc_64 <- c()
for (i in 1:n_64) { 
  c_64 <- rbind(c_64,rawDat_64$Contribution[rawDat_64$uniqueSubno==i])
  Ga_64 <- rbind(Ga_64,rawDat_64$averageC[rawDat_64$uniqueSubno==i])
  Gc_64 <- rbind(Gc_64,rawDat_64$sumOthersContrib[rawDat_64$uniqueSubno==i])
}

# compile data from each condition into 4D matrix
c <- array(0,c(n_16,ntrials,2))
c[,,1] <- c_16
c[,,2] <- c_64

Ga <- array(0,c(n_16,ntrials,2))
Ga[,,1] <- Ga_16
Ga[,,2] <- Ga_64

Gc <- array(0,c(n_16,ntrials,2))
Gc[,,1] <- Gc_16
Gc[,,2] <- Gc_64

c_choice_index <- c

data <- list("groupSize","n_16", "ntrials", "ntokens", "pi", "vals","c","Gc","c_choice_index","Ga") #data inputted into jags
params <- c("mu_alpha","mu","alpha","c","sigma") #parameters we'll track in jags

samples <- jags(data, inits=NULL, params,
                model.file ="multiplicationFactor_ttest.txt",
                n.chains=3, n.iter=15000, n.burnin=1000, n.thin=3)


# calculate Bayes factors for difference using logspline fit
prior <- dnorm(0,1)
fit.posterior <- logspline(samples$BUGSoutput$sims.list$mu_alpha)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior/posterior

# credible intervals
qlogspline(0.025,fit.posterior) #2.5% CI
qlogspline(0.975,fit.posterior) #97.5% CI

# convert posteriors back into probability space
pCC_64 <- pnorm(samples$BUGSoutput$sims.list$mu+samples$BUGSoutput$sims.list$alpha)
pCC_16 <- pnorm(samples$BUGSoutput$sims.list$mu)

# find the MAPs for each individual subject
MAP_16 <- array(0,c(n_16))
MAP_64 <- array(0,c(n_16))
CI_16 <- array(0,c(n_16,2))
CI_64 <- array(0,c(n_16,2))

for (i in 1:n_16) {
  MAP_16[i] <- density(pCC_16[,i])$x[which(density(pCC_16[,i])$y==max(density(pCC_16[,i])$y))]
  MAP_64[i] <- density(pCC_64[,i])$x[which(density(pCC_64[,i])$y==max(density(pCC_64[,i])$y))]
  
}

########################################################################

par(mfrow=c(1,3),cex=1.5)

#make a plot to show trends in raw contributions
plot(colMeans(c_16),ylim=c(0,40),xlim=c(0,20),type='l',col="red",lty=1,lwd=2,
     xlab = "Trial", ylab = "Investment amount", main = expression("Effects of Payoff"),
     axes = FALSE)
axis(1)
axis(2)
# standard deviations above
lines(colMeans(c_16)+(colMeans(abs(c_16-colMeans(c_16))))/sqrt(n_16),lty=3,col="red",lwd=2)
# standard deviations below
lines(colMeans(c_16)-(colMeans(abs(c_16-colMeans(c_16))))/sqrt(n_16),lty=3,col="red",lwd=2)

# mean for 64 condition
lines(colMeans(c_64),lty=1,col="blue",lwd=2)
# standard deviations above
lines(colMeans(c_64)+(colMeans(abs(c_64-colMeans(c_64))))/sqrt(n_64),lty=3,col="blue",lwd=2)
# standard deviations below
lines(colMeans(c_64)-(colMeans(abs(c_64-colMeans(c_64))))/sqrt(n_64),lty=3,col="blue",lwd=2)

legend (x = 7, y = 20, legend = c(expression(paste(pi,"= 6.4")), expression(paste(pi,"= 1.6"))), 
        col = c("blue","red"), bty = "n", lwd = 2,cex=.8)


# make a plot to posterior for difference
plot(density(samples$BUGSoutput$sims.list$mu_alpha),xlab = expression(paste(alpha[mu])),
     ylab = "",lwd = 3, frame=FALSE,main=expression(paste("Posterior for ",alpha[mu]))) 


Condition <- c(rep(1,n_16),rep(2,n_16))
MAPs <- c(MAP_16,MAP_64)

plot(jitter(Condition,.2),MAPs,xlim=c(0.5,2.5),ylim=c(0,1),ylab = "",
  frame=FALSE,main=expression("Prob. Coop."),axes=FALSE,xlab=expression(pi))
  axis(1, at=1:2, labels=c(1.6,6.4))
  axis(2)

for (i in 1:n_16) {
  lines(c(MAP_16[i],MAP_64[i]),lty=3)
}