
# Example code from Kery & Schaub, implementing Conn & Cooch 2009

library(jagsUI)
library(R2OpenBUGS)

# Function to simulate multistate capture-recapture data
simul.me <- function(PSI.STATE, PSI.OBS, PSI.INI, OBS.INI, marked, unobservable = NA){
   # Unobservable: number of state that is unobservable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked)
   
   for (i in 1:sum(marked)){
      # Initial state
      CH.TRUE[i,mark.occ[i]] <- which(rmultinom(1, 1, PSI.INI[,i,mark.occ[i]])==1)
      CH[i,mark.occ[i]] <- which(rmultinom(1, 1, OBS.INI[CH.TRUE[i,mark.occ[i]],,i,mark.occ[i]])==1)
      for (t in (mark.occ[i]+1):n.occasions){
         # Multinomial trials for state transitions
         if (mark.occ[i]==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Multinomial trials for observation process
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Replace the NA and the highest state number (dead) in the file by 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.OBS)[2]] <- 0
   CH[CH==unobservable] <- 0
   id <- numeric(0)
   for (i in 1:dim(CH)[1]){
      z <- min(which(CH[i,]!=0))
      ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      }
   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
   # CH: capture histories to be used
   # CH.TRUE: capture histories with perfect observation
   }


# # Example 1. A model where the true state of an individual may not always be seen
# 
# 
# # Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
# phiA <- 0.8
# phiB <- 0.6
# psiAB <- 0.3
# psiBA <- 0.5
# pA <- 0.7
# pB <- 0.4
# dA <- 0.9     # probability that the state of an individual that is in state A is correctly observed 
# dB <- 0.6     # probability that the state of an individual that is in state B is correctly observed
# psi <- 0.75   # probability that an individual is in state A when marked
# 
# n.occasions <- 6
# n.states <- 3
# n.obs <- 4
# marked <- rep(2500, n.occasions)  
# 
# # Define matrices with survival, transition and recapture probabilities
# # These are 4-dimensional matrices, with 
#    # Dimension 1: state of departure
#    # Dimension 2: state of arrival
#    # Dimension 3: individual
#    # Dimension 4: time
# 
# # 1. State process matrix
# totrel <- sum(marked)
# PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
# for (i in 1:totrel){
#    for (t in 1:(n.occasions-1)){
#       PSI.STATE[,,i,t] <- matrix(c(
#       phiA*(1-psiAB), phiA*psiAB,     1-phiA,
#       phiB*psiBA,     phiB*(1-psiBA), 1-phiB,
#       0,              0,              1       ), nrow = n.states, byrow = TRUE)
#       } #t
#    } #i
# 
# # 2.Observation process matrix
# PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
# for (i in 1:totrel){
#    for (t in 1:(n.occasions-1)){
#       PSI.OBS[,,i,t] <- matrix(c(
#       dA * pA, 0,  (1-dA) * pA, 1-pA,
#       0,  dB * pB, (1-dB) * pB, 1-pB,
#       0,  0,  0, 1       ), nrow = n.states, byrow = TRUE)
#       } #t
#    } #i
# 
# # 3. Matrix with the initial state probabilities
# PSI.INI <- array(NA, dim=c(n.states-1, totrel, n.occasions))
# for (i in 1:totrel){
#    for (t in 1:n.occasions){
#       PSI.INI[1,i,t] <- psi
#       PSI.INI[2,i,t] <- 1-psi
#       }
#    }
# 
# # 4. Assignment matrix at initial capture
# OBS.INI <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions))
# for (i in 1:totrel){
#    for (t in 1:n.occasions){
#       OBS.INI[,,i,t] <- matrix(c(
#       dA , 0,  1-dA, 0,
#       0,  dB,  1-dB, 0,
#       0,  0,  0, 1       ), nrow = n.states, byrow = TRUE)
#       } #t
#    } #i
# 
# # Execute function
# sim <- simul.me(PSI.STATE, PSI.OBS, PSI.INI, OBS.INI, marked)
# CH <- sim$CH
# table(CH) # 0 = not seen, 1 = seen alive in A, 2 = seen alive in B, 3 = seen alive but unknown site ("state") 
# 
# 
# # Compute vector with occasion of first capture
# get.first <- function(x) min(which(x!=0))
# f <- apply(CH, 1, get.first)
# 
# # Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# # 1 = seen alive in A, 2 = seen alive in B, 3 = seen alive, unknown site, 4 = not seen
# rCH <- CH          # Recoded CH
# rCH[rCH==0] <- 4
# 
# 
# # Specify model in jags
# 
# model <- function () {
# 
# # -------------------------------------------------
# # Parameters:
# # phiA: survival probability at site A
# # phiB: survival probability at site B
# # psiAB: movement probability from site A to site B
# # psiBA: movement probability from site B to site A
# # pA: recapture probability at site A
# # pB: recapture probability at site B
# # dA: probability to correctly assign site A to an observed individual that is at site A
# # dB: probability to correctly assign site B to an observed individual that is at site B
# # pi: probability that an individual is at site A when marked 
# # -------------------------------------------------
# # States (S):
# # 1 alive at A
# # 2 alive at B
# # 3 dead
# # Observations (O):  
# # 1 seen at A 
# # 2 seen at B
# # 3 seen, but site unknown
# # 4 not seen
# # -------------------------------------------------
# 
# # Priors and constraints
# for (i in 1:nind){
#    for (t in 1:(n.occasions-1)){
#       phiA[i,t] <- mean.phi[1]
#       phiB[i,t] <- mean.phi[2]
#       psiAB[i,t] <- mean.psi[1]
#       psiBA[i,t] <- mean.psi[2]
#       pA[i,t] <- mean.p[1]
#       pB[i,t] <- mean.p[2]
#       dA[i,t] <- mean.d[1]
#       dB[i,t] <- mean.d[2]
#       pi[i,t] <- mean.pi
#       } # t
#    } # i
# 
# for (u in 1:2){
#    mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
#    mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
#    mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
#    mean.d[u] ~ dunif(0, 1)      # Priors for mean assignment error
#    }
# mean.pi ~ dunif(0, 1)
# 
# # Define state-transition and observation matrices
# for (i in 1:nind){
#    for (t in f[i]:(n.occasions-1)){
#       # Define state probability when an individual is marked
#       pj[i,t,1] <- pi[i,t]
#       pj[i,t,2] <- 1-pi[i,t]
#       
#       # Define probabilities of O(t) given S(t) at the occasion of marking (this corresponds to B(0) in E-SURGE notation)
#       po1[1,i,t,1] <- dA[i,t]
#       po1[1,i,t,2] <- 0
#       po1[1,i,t,3] <- 1-dA[i,t]
#       po1[1,i,t,4] <- 0
#       po1[2,i,t,1] <- 0
#       po1[2,i,t,2] <- dB[i,t]
#       po1[2,i,t,3] <- 1-dB[i,t]
#       po1[2,i,t,4] <- 0
#       po1[3,i,t,1] <- 0
#       po1[3,i,t,2] <- 0
#       po1[3,i,t,3] <- 0
#       po1[3,i,t,4] <- 1
#    
#       # Define probabilities of state S(t+1) given S(t)   
#       ps[1,i,t,1] <- phiA[i,t] * (1-psiAB[i,t])
#       ps[1,i,t,2] <- phiA[i,t] * psiAB[i,t]
#       ps[1,i,t,3] <- 1-phiA[i,t]
#       ps[2,i,t,1] <- phiB[i,t] * psiBA[i,t]
#       ps[2,i,t,2] <- phiB[i,t] * (1-psiBA[i,t])
#       ps[2,i,t,3] <- 1-phiB[i,t]
#       ps[3,i,t,1] <- 0
#       ps[3,i,t,2] <- 0
#       ps[3,i,t,3] <- 1
#       
#       # Define probabilities of O(t) given S(t)
#       po[1,i,t,1] <- dA[i,t] * pA[i,t]
#       po[1,i,t,2] <- 0
#       po[1,i,t,3] <- (1-dA[i,t]) * pA[i,t]
#       po[1,i,t,4] <- 1-pA[i,t]
#       po[2,i,t,1] <- 0
#       po[2,i,t,2] <- dB[i,t] * pB[i,t]
#       po[2,i,t,3] <- (1-dB[i,t]) * pB[i,t]
#       po[2,i,t,4] <- 1-pB[i,t]
#       po[3,i,t,1] <- 0
#       po[3,i,t,2] <- 0
#       po[3,i,t,3] <- 0
#       po[3,i,t,4] <- 1
#       } #t
#    } #i
# 
# # Likelihood 
# for (i in 1:nind){
#    # Define latent state at first capture
#    z[i,f[i]] ~ dcat(pj[i, f[i],])
#    y[i,f[i]] ~ dcat(po1[z[i,f[i]], i, f[i],])
#    for (t in (f[i]+1):n.occasions){
#       # State process: draw S(t) given S(t-1)
#       z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
#       # Observation process: draw O(t) given S(t)
#       y[i,t] ~ dcat(po[z[i,t], i, t-1,])
#       } #t
#    } #i
# }
# 
# # Function to create initial values for unknown z
# me.init.z <- function(ch, f){
#    g <- which(ch == 3)
#    g <- c(g, which(ch==4))
#    states <- max(ch, na.rm = TRUE)
#    known.states <- 1:(states-2)
#    ch[g] <- sample(known.states, length(g), replace = TRUE)   
#    for (i in 1:nrow(ch)){
#       fi <- ch[i,f[i]]
#       ch[i,1:f[i]] <- NA
#       ch[i,f[i]] <- fi
#       }
#    return(ch)
#    }
# write.model(model, "model.txt")
# model.file = paste(getwd(),"model.txt", sep="/")
# 
# 
# # Bundle data
# jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(CH)[1])
# 
# # Initial values
# inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), 
#                          mean.p = runif(2, 0, 1), z = me.init.z(rCH, f))}  
# 
# # Parameters monitored
# parameters <- c("mean.phi", "mean.psi", "mean.p", "mean.d", "mean.pi")
# 
# # MCMC settings
# nc <- 3; nAdapt <- 1000; nb <- 2000; ni <- 5000+nb; nt <- 1  
# 
# # Call jagsUI
# me <- jags(jags.data, inits, parameters, n.chains = nc, model.file = model.file,
#            n.adapt=nAdapt, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=F)
# 
# print(me, digits = 3)
# 

##############################################
#
# Example 2. A model where the disease may not always be seen
#

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiH <- 0.8  # survial if healthy
phiD <- 0.6  # survival if diseased
psiHD <- 0.3 # probability healthy individual transitions to disease state, conditional on survival
psiDH <- 0.5 # probability diseased individual transitions to healthy state, conditional on survival
pH <- 0.7    # probability of detection if in healthy state
pD <- 0.4    # probability of detection if in diseased state
b <- 0.2     # probability that the disease is not detected
psi <- 0.75  # probability that an individual is in healthy state when marked

n.occasions <- 4
n.states <- 3
n.obs <- 3
marked <- rep(10, n.occasions)  

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phiH*(1-psiHD), phiH*psiHD,     1-phiH,
      phiD*psiDH,     phiD*(1-psiDH), 1-phiD,
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      pH, 0, 1-pH,
      b * pD,  (1-b) * pD, 1-pD,
      0,  0, 1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 3. Matrix with the initial state probabilities
PSI.INI <- array(NA, dim=c(n.states-1, totrel, n.occasions))
for (i in 1:totrel){
   for (t in 1:n.occasions){
      PSI.INI[1,i,t] <- psi
      PSI.INI[2,i,t] <- 1-psi
      }
   }

# 4. Assignment matrix at initial capture
OBS.INI <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions))
for (i in 1:totrel){
   for (t in 1:n.occasions){
      OBS.INI[,,i,t] <- matrix(c(
      1, 0,  0,
      b,  1-b,  0,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i


# Execute function
sim <- simul.me(PSI.STATE, PSI.OBS, PSI.INI, OBS.INI, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen without disease, 2 = seen with disease, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3


# Specify model in BUGS language

model_conn <- function () {

# -------------------------------------------------
# Parameters:
# phiH: survival probability when not infected by disease
# phiD: survival probability when infected by disease
# psiHD: probability that a healty individual gets the disease 
# psiDH: probability that a ill individual gets healthy again
# pH: recapture probability of healthy individuals
# pD: recapture probability of infected individuals
# b: probability that the disease is not detected
# pi: probability that an individual is healthy when marked (first encounter) 
# -------------------------------------------------
# States (S):
# 1 alive not infected
# 2 alive with disease
# 3 dead
# Observations (O):  
# 1 seen without disease 
# 2 seen with disease
# 3 not seen
# -------------------------------------------------

# Priors and constraints
for (i in 1:nind){
   for (t in 1:(n.occasions-1)){
      phiH[i,t] <- mean.phi[1]
      phiD[i,t] <- mean.phi[2]
      psiHD[i,t] <- mean.psi[1]
      psiDH[i,t] <- mean.psi[2]
      pH[i,t] <- mean.p[1]
      pD[i,t] <- mean.p[2]
      b[i,t] <- mean.b
      pi[i,t] <- mean.pi
      } # t
   } # i

for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }
mean.b ~ dunif(0, 1)      # Priors for mean assignment error
mean.pi ~ dunif(0, 1)

# Define state-transition and observation matrices
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      # Define state probability when an individual is marked
      pj[i,t,1] <- pi[i,t]
      pj[i,t,2] <- 1-pi[i,t]
      
      # Define probabilities of O(t) given S(t) at the occasion of marking (this corresponds to B(0) in E-SURGE notation)
      po1[1,i,t,1] <- 1
      po1[1,i,t,2] <- 0
      po1[1,i,t,3] <- 0
      po1[2,i,t,1] <- b[i,t]
      po1[2,i,t,2] <- 1-b[i,t]
      po1[2,i,t,3] <- 0
      po1[3,i,t,1] <- 0
      po1[3,i,t,2] <- 0
      po1[3,i,t,3] <- 1
   
      # Define probabilities of state S(t+1) given S(t)   
      ps[1,i,t,1] <- phiH[i,t] * (1-psiHD[i,t])
      ps[1,i,t,2] <- phiH[i,t] * psiHD[i,t]
      ps[1,i,t,3] <- 1-phiH[i,t]
      ps[2,i,t,1] <- phiD[i,t] * psiDH[i,t]
      ps[2,i,t,2] <- phiD[i,t] * (1-psiDH[i,t])
      ps[2,i,t,3] <- 1-phiD[i,t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pH[i,t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pH[i,t]
      po[2,i,t,1] <- b[i,t] * pD[i,t]
      po[2,i,t,2] <- (1-b[i,t]) * pD[i,t]
      po[2,i,t,3] <- 1-pD[i,t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] ~ dcat(pj[i, f[i],])
   y[i,f[i]] ~ dcat(po1[z[i,f[i]], i, f[i], ])
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}

write.model(model_conn, "model.txt")
model.file = paste(getwd(),"model.txt", sep="/")

# Function to create initial values for unknown z
me.init.z <- function(ch, f){
   g <- which(ch == 3)
   g <- c(g, which(ch==4))
   states <- max(ch, na.rm = TRUE)
   known.states <- 1:(states-2)
   ch[g] <- sample(known.states, length(g), replace = TRUE)   
   for (i in 1:nrow(ch)){
      fi <- ch[i,f[i]]
      ch[i,1:f[i]] <- NA
      ch[i,f[i]] <- fi
      }
   return(ch)
   }


# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), 
                         mean.psi = runif(2, 0, 1), 
                         mean.p = runif(2, 0, 1), z = me.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p", "mean.b", "mean.pi")

# MCMC settings
nc <- 3; nAdapt <- 100; nb <- 200; ni <- 500+nb; nt <- 1  

# Call jagsUI
me <- jags(jags.data, inits, parameters, model.file = model.file, n.chains = nc, 
           n.adapt=nAdapt, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

print(me, digits = 3)
