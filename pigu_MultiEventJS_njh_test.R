
# Following example code from Kery & Schaub, implementing Conn & Cooch 2009

rm(list=ls())
library(jagsUI)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("~/Documents/SAFS/PigeonGuillemots")
#setwd("C:\\Nathan\\UW\\projects\\PigeonGuillemotNestSurvival\\data")


#load data from PigeonGuillemot_whidbey.Rmd

burrow_CH <- read.csv("burrow_CH.csv", header = T, stringsAsFactors = F)


data_MEJS <- burrow_CH %>%
  select(region, year, site, study_day, burrow_name, capt_hist) %>%
  distinct()
 


####-#### reduce data for testing

data_MEJS <-subset(data_MEJS, site == "Double Bluff North" &  year == 2009)

test <- data_MEJS[duplicated(data_MEJS),]


data_MEJS_wide <- data_MEJS %>% 
  distinct(region, year, site, burrow_name, study_day, capt_hist) %>% 
  dcast(region + year + site + burrow_name ~ study_day, value.var = 'capt_hist', fill = 4, fun.aggregate = mean)

CH_MEJS <- as.matrix(data_MEJS_wide %>% select(-c(region, year, site, burrow_name)), rownames.force = F)
colnames(CH_MEJS) <- NULL
CH_MEJS <- round(CH_MEJS, 0)

## Create binary effort matrix with 1 if site surveyed on day t, 0 otherwise.
EFFORT_MEJS <- matrix(0, nr=dim(CH_MEJS)[1], nc=dim(CH_MEJS)[2])
EFFORT_MEJS[!is.na(CH_MEJS)] <- 1

## check out some burrow specific data
junk<-matrix(0, nr=3, nc=dim(CH_MEJS)[2])
for(j in 1:dim(junk)[2]){
  junk[as.numeric(names(table(CH_MEJS[,j]))),j]<-table(CH_MEJS[,j])
}

# Specify model in jags

sink("model_MEJS.txt")
cat("
model {
# -------------------------------------------------
# Parameters:
# phiA: survival probability from egg to chick
# phiB: survival probability from chick to fledge
# psiAB: probability of transitioning from egg to chick
# pA: detection probability of egg burrow
# pB: detection probability of chick burrow
# b: conditional on observing a visit to a chick burrow, probability the delivery was a 'burrow visit' 
# gamma: entry probability
# alpha: conditional on entry at occasion 1, probability burrow had a chick. alpha set to 0 for t>1. So if a burrow is initiated after day one, it must start as an egg burrow.  

# -------------------------------------------------
# States (S):
# 1 not entered
# 2 alive as egg
# 3 alive as chick
# 4 terminated; dead or fledged
# Observations (O):  
# 1 Burrow visit
# 2 Prey visit
# 3 not seen
# -------------------------------------------------

# Priors and constraints
b~dunif(0,1)            ## start with simple constant model
mean.psiAB~dunif(0,1)   ## start with simple constant model
alpha[1] ~ dunif(0,1)   ## if entery at occasion1, probability entered as chick burrow
alphaKeep <-alpha[1]    ## this is the only one we want to keep 

for (u in 1:2){
 mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
 mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
}
for (t in 1:(n.occasions-1)){
 phiA[t]     <- mean.phi[1]   # egg survival
 phiB[t]     <- mean.phi[2]   # chick survival
 gamma[t]    ~ dunif(0, 1)    # Prior for entry probabilities at occasion t
 pA[t]       <- mean.p[1]     # egg burrow detection
 pB[t]       <- mean.p[2]     # chick burrow detection
 psiAB[t]    <- mean.psiAB
 alpha[t+1]  <-0              #probability of entry as chick on occasion t+1 
}

# Likelihood 
for (i in 1:M){
 # Define latent state at first occasion
 z[i,1] <- 1            # Make sure that all M individuals are in state 1 at t=1 (e.g., one day prior to season)
 for (t in 2:n.occasions){
  # State process: draw S(t) given S(t-1)
  z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
  # Observation process: draw O(t) given S(t)
  y[i,t] ~ dcat(po[z[i,t], i, t-1,])  
   #this follow Kery and Schaub 10.3.2, indexing is tricky. e.g., y[i,2] is a function of state[i,2] and observation array t=1
   #so if detection covariates are used, DO NOT ADD DUMMY OCCASION to detection covariates as detection data at occasion t uses array t-1

  # Derive some stuff along the way
  egg[i,t-1]    <- equals(z[i,t], 2)
  chick[i,t-1]  <- equals(z[i,t], 3)
  active[i,t-1] <- max(egg[i,t-1],chick[i,t-1]) #egg or chick
 } #t
 everEgg[i]    <- max(egg[i,])
 everChick[i]  <- max(chick[i,])
 EverActive[i] <- max(everChick[i],everEgg[i])
} #i


# Define transition and observation matrices
 for (i in 1:M){
  for (t in 1:(n.occasions-1)) {
      # Define probabilities of state S(t+1) given S(t)   
      ps[1,i,t,1] <- 1-gamma[t]               #probability of not entering
      ps[1,i,t,2] <- gamma[t]*(1-alpha[t])    #probability of entering as egg
      ps[1,i,t,3] <- gamma[t]*alpha[t]        #probability of entering as chick. Note: only alpha[1]>0
      ps[1,i,t,4] <- 0                        #probability of entering as terminated
      ps[2,i,t,1] <- 0                        #probability egg goes to 'not entered'
      ps[2,i,t,2] <- phiA[t]*(1-psiAB[t])     #probability of surviving egg state and not transitioning
      ps[2,i,t,3] <- phiA[t]*psiAB[t]         #probability of surviving egg state and hatching
      ps[2,i,t,4] <- 1-phiA[t]                #probability of a failed egg
      ps[3,i,t,1] <- 0                        #probability chick goes to 'not entered'
      ps[3,i,t,2] <- 0                        #probability chickgoes to egg
      ps[3,i,t,3] <- phiB[t]                  #probability of surviving chick stage
      ps[3,i,t,4] <- 1-phiB[t]                #probability of failed chick
      ps[4,i,t,1] <- 0                        #probability terminated goes to 'not entered'
      ps[4,i,t,2] <- 0                        #probability terminated goes to egg (maybe?)
      ps[4,i,t,3] <- 0                        #probability terminated goes to chick
      ps[4,i,t,4] <- 1                        #probability terminated goes to terminated
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0                        #probability an 'not entered' burrow is detected with a burrow visit
      po[1,i,t,2] <- 0                        #probability an 'not entered' burrow is detected with a prey visit
      po[1,i,t,3] <- 1                        #probability an 'not entered' burrow is not detected
      po[2,i,t,1] <- pA[t]*eff[i,t]           #probability an egg burrow is detected with a burrow visit
      po[2,i,t,2] <- 0                        #probability an egg burrow is detected with a prey visit
      po[2,i,t,3] <- 1-pA[t]*eff[i,t]         #probability an egg burrow is not detected
      po[3,i,t,1] <- b * pB[t]*eff[i,t]       #probability a chick burrow is detected with a burrow visit 
      po[3,i,t,2] <- (1 - b) * pB[t]*eff[i,t] #probability a chick burrow is detected with a prey visit
      po[3,i,t,3] <- 1 - pB[t]*eff[i,t]       #probability a chick burrow is not detected
      po[4,i,t,1] <- 0                        #probability a terminated burrow is detected with a burrow visit
      po[4,i,t,2] <- 0                        #probability a terminated burrow is detected with a prey visit
      po[4,i,t,3] <- 1                        #probability a terminated burrow is not detected
  } #t
 }#M
#derive abundances (lots of options here, e.g., site specific, annual, etc.) 
 for (t in 1:(n.occasions-1)){
  N.egg[t]    <- sum(egg[1:M,t])
  N.chick[t]  <- sum(chick[1:M,t])
  N.active[t] <- sum(active[1:M,t])
  qgamma[t] <- 1-gamma[t]
  birthProb[t] <- cprob[t] / psi      # Entry probability
}
 cprob[1] <- gamma[1]
 for (t in 2:(n.occasions-1)){
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
 }#t
 psi <- sum(cprob[])            # Inclusion probability

Nstar<-sum(EverActive[1:M])
 
}#model
",fill=TRUE)
sink()



##########################
##-## Data augmentation may need to be site x year specific. This will take some thought and care in developing augmented datasets
   
# Augment data
nz <- 100
CH.aug <- rbind(CH_MEJS, matrix(3, nrow = nz, ncol = dim(CH_MEJS)[2])) #3 for not seen

#augmented effort (This will take some thought...)
eff.aug<- rbind(EFFORT_MEJS, matrix(EFFORT_MEJS[1,], nrow=nz, ncol=dim(EFFORT_MEJS)[2], byrow=T))
eff.aug<-eff.aug[,-1] #remove first column (dummy column). See note in JAGS model
##########################

     
# Bundle data
jags.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1], eff=eff.aug)

get.first <- function(x) min(which(x < 3))
f <- apply(CH_MEJS, 1, get.first)
    
get.last <- function(x) max(which(x < 3))
l <- apply(CH_MEJS, 1, get.last)

#inits <- CH_MEJS
js.me.init <- function(ch, f, l, nz){
  inits <- matrix(NA, nr=dim(ch)[1], nc=dim(ch)[2])
  inits[,1] <- NA  #deterministic on occasion 1

  for(i in 1:dim(ch)[1]) { 
    inits[i,1:f[i]-1] <- 1  #1 before first capture 
    if(l[i]<dim(ch)[2]){
     inits[i, (l[i]+1):dim(ch)[2]] <- rep(3, (dim(ch)[2]-l[i])) #3 from last obs to end
    }
    if(sum(ch[i,] == 2, na.rm=T)>0){ #sites with observed prey delivery
       if(f[i]==min(which(ch[i,]  == 2))){  #if first obs is a prey delivery
        inits[i, f[i]:l[i]]<-3              #chick until end of study 
        if(f[i]>2) inits[i, 2:(f[i]-1)]<-2  #egg prior to chick

       }else{
        inits[i, f[i]:(min(which(ch[i,]  == 2))-1)]<-2
        inits[i, min(which(ch[i,]  == 2)):l[i]]<-3
       }
    }else{ #never detected a prey delivery
      if(sum(ch[i,] == 2, na.rm=T)==0){      
        inits[i, f[i]:(max(which(ch[i,]  == 1)))]<-2
      }
    }
   }#i
  return(inits)
}


zInitsObs = js.me.init(CH_MEJS, f, l, nz)
zInits<-rbind(zInitsObs,matrix(1,nr=nz, nc=dim(CH_MEJS)[2]))
zInits[,1]<-NA


inits <- function(){list(mean.phi = runif(2, 0, 1), 
                         mean.p = runif(2, 0, 1), 
                         z = zInits)}    
     
# Parameters monitored
parameters <- c("mean.phi", "mean.p", "b", "gamma", "alphaKeep", "mean.psiAB","N.egg","N.chick", "N.active","Nstar","psi","birthProb")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 1000
nc <- 3
     
js_me <- jags(jags.data, inits, parameters, model.file = "model_MEJS.txt",
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
     
print(js_me, digits = 3)

plot(js_me$samples[,"mean.phi[2]"])

outmat<-as.matrix(js_me$samples)
plot(apply(outmat[,grep("N.active\\[", colnames(outmat))],2,median), pch=16, type="b", ylab="Number", xlab=c("Day of season"))
points(apply(outmat[,grep("N.egg\\[", colnames(outmat))],2,median), pch=16, type="b", col="red")
points(apply(outmat[,grep("N.chick\\[", colnames(outmat))],2,median), pch=16, type="b", col="blue")
legend("topleft", legend=c("Active","Egg", "Chick"), text.col=c("black", "red", "blue"), bty="n")

