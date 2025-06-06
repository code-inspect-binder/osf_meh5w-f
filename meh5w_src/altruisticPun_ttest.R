set.seed(1982)
library(R2jags)
library(polspline)

# This script applies a hierarchical t-test over latent mixture model parameters
# to data from the paper Hunger affects social decisions...
# Pre-print available https://psyarxiv.com/67abq/
# Data available 
# https://figshare.com/articles/Data_and_code_for_Fraser_and_Nettle_Hunger_affects_social_decisions_in_a_Public_Goods_Game_but_not_an_Ultimatum_Game_/11807346/1

# analysis is not conducted on hunger condition. Data are used to analyse only effects
# of altruistic punishment

nsubs <- 4
groupSize <- nsubs
ntrials <- 10
pi <- 1.4
ntokens <- 20
vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

setwd("SET DATA DIRECTORY")

rawDat <- read.csv("HungerAffectsSocialDecisions.csv") # Public goods game

group_names <- unique(rawDat$UniqueGroup)
ngroups <- length(group_names)

# create data matrices - subject x trial x group x condition

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

# data for punishment condition #
c_punish <- array(0,c(nsubs,ntrials,ngroups))
Ga_punish <- array(0,c(ntrials,ngroups))
Gc_punish <- array(0,c(nsubs,ntrials,ngroups))

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="Y"][1:10],
                    rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="Y"][11:20],
                    rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="Y"][21:30],
                    rawDat$Contribution[rawDat$UniqueGroup==group_names[g]&rawDat$Punishment.Round=="Y"][31:40])
  
  Ga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:nsubs) {
    Gc_punish[,,g] <- colSums(c_punish[-s,,g])
  }
}

# compile data from each condition into 4D matrix
c <- array(0,c(nsubs,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Ga <- array(0,c(ntrials,ngroups,2))
Ga[,,1] <- Ga_no_punish
Ga[,,2] <- Ga_punish

Gc <- array(0,c(nsubs,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

c_choice_index <- c

data <- list("groupSize", "ngroups", "ntrials", "ntokens", "pi", "vals","c","Gc","c_choice_index","Ga") #data inputted into jags
params <- c("alpha","mu","mu_alpha") #parameters we'll track in jags

samples <- jags(data, inits=NULL, params,
                model.file ="altruisticPun_ttest.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1)

# visual difference
plot(density(samples$BUGSoutput$sims.list$mu_alpha),xlab = expression(paste(alpha)),lwd = 3, main = "",frame=FALSE) 

# calculate Bayes factors for difference using logspline fit
prior <- dnorm(0,1)
fit.posterior <- logspline(samples$BUGSoutput$sims.list$mu_alpha)
posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior/posterior

# credible intervals
qlogspline(0.025,fit.posterior) #2.5% CI
qlogspline(0.975,fit.posterior) #97.5% CI

# convert posteriors back into probability space
pCC_punish <- pnorm(samples$BUGSoutput$sims.list$mu+samples$BUGSoutput$sims.list$alpha)
pCC_noPunish <- pnorm(samples$BUGSoutput$sims.list$mu)

# find the MAPs for each individual subject
MAP_noPunish <- array(0,c(nsubs,ngroups))
MAP_punish <- array(0,c(nsubs,ngroups))
CI_noPunish <- array(0,c(n_16,2))
CI_punish <- array(0,c(n_16,2))

for (g in 1:ngroups) {
  for (s in 1:nsubs) {
    MAP_punish[s,g] <- density(pCC_punish[,s,g])$x[which(density(pCC_punish[,s,g])$y==max(density(pCC_punish[,s,g])$y))]
    MAP_noPunish[s,g] <- density(pCC_noPunish[,s,g])$x[which(density(pCC_noPunish[,s,g])$y==max(density(pCC_noPunish[,s,g])$y))]
  
  }
}

MAP_noPunish <- c(MAP_noPunish[1,],MAP_noPunish[2,],MAP_noPunish[3,],MAP_noPunish[4,])
MAP_punish <- c(MAP_punish[1,],MAP_punish[2,],MAP_punish[3,],MAP_punish[4,])

########################################################################

#make arrays of individual subs for plotting
c_punish_nogroups <- t(cbind(c_punish[1,,],c_punish[2,,],c_punish[3,,],c_punish[4,,]))
c_nopunish_nogroups <- t(cbind(c_no_punish[1,,],c_no_punish[2,,],c_no_punish[3,,],c_no_punish[4,,]))

par(mfrow=c(1,3),cex=1.5)

#make a plot to show trends in raw contributions
plot(colMeans(c_punish_nogroups),ylim=c(0,20),xlim=c(0,10),type='l',col="red",lty=1,lwd=2,
     xlab = "Trial", ylab = "Investment amount", main = expression("Effects of Punishment"),
     axes = FALSE)
axis(1)
axis(2)
# standard deviations above
lines(colMeans(c_punish_nogroups)+(colMeans(abs(c_punish_nogroups-colMeans(c_punish_nogroups))))/sqrt(nsubs*ngroups),lty=3,col="red",lwd=2)
# standard deviations below
lines(colMeans(c_punish_nogroups)-(colMeans(abs(c_punish_nogroups-colMeans(c_punish_nogroups))))/sqrt(nsubs*ngroups),lty=3,col="red",lwd=2)

# mean for 64 condition
lines(colMeans(c_nopunish_nogroups),lty=1,col="blue",lwd=2)
# standard deviations above
lines(colMeans(c_nopunish_nogroups)+(colMeans(abs(c_nopunish_nogroups-colMeans(c_nopunish_nogroups))))/sqrt(nsubs*ngroups),lty=3,col="blue",lwd=2)
# standard deviations below
lines(colMeans(c_nopunish_nogroups)-(colMeans(abs(c_nopunish_nogroups-colMeans(c_nopunish_nogroups))))/sqrt(nsubs*ngroups),lty=3,col="blue",lwd=2)

legend (x = 1, y = 20, legend = c("Punish", "Standard"), 
        col = c("red","blue"), bty = "n", lwd = 2,cex =.8)


# make a plot to posterior for difference
plot(density(samples$BUGSoutput$sims.list$mu_alpha),xlab = expression(paste(alpha[mu])),ylab = "",
     lwd = 3, frame=FALSE,main=expression(paste("Posterior for ",alpha[mu]))) 


Condition <- c(rep(1,ngroups*nsubs),rep(2,ngroups*nsubs))
MAPs <- c(MAP_noPunish,MAP_punish)

plot(jitter(Condition,.2),MAPs,xlim=c(0.5,2.5),ylim=c(0,1),
     frame=FALSE,main=expression("Prob. Coop."),axes=FALSE,xlab="Punishment",ylab = "")
axis(1, at=1:2, labels=c("No","Yes"))
axis(2)

for (i in 1:ngroups*nsubs) {
  lines(c(MAP_noPunish[i],MAP_punish[i]),lty=3)
}

