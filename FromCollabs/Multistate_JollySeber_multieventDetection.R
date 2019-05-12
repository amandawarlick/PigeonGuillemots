rm(list=ls())
library(jagsUI)

# Define parameter values
n.occasions <- 50       # Number of capture occasions
Nstar <- 1000           # Superpopulation size (total number of burrows active on at least one day)
n.states <- 4           # not entered, egg, chick, terminated
n.obs <- 3              # burrow visit, prey delivery, not detected
b <- c(.5, seq(.3, .05, length = n.occasions-1))      # Entry probabilities decrease with time
alpha <- c(.40, rep(0,n.occasions-1))    # Prob burrow is a "chick burrow" given entry (only occurs on day 1)
phiA <- 0.8             # survival state A (egg burrow)
phiB <- 0.5             # survival state B (chick burrow)
psiAB <- 0.25           # transition A to B conditional on survival
pA <- 0.5               # detection for egg burrow
pB <- 0.75              # detection for chick burrow
beta <- 0.25            # conditional on observing a visit to a chick burrow, probability the delivery was a 'burrow visit'

#Define matrics

#state at entry
INIT.STATE <- matrix(NA, nr = n.occasions, nc = n.states) 
 for(t in 1:n.occasions){
  INIT.STATE[t,] <- c(0, 1-alpha[t], alpha[t], 0) #conditional on entry at time t (b), prob of entry as chick burrow 
 }

#state process, conditional on entry 
STATE_MATRIX <- array(NA, dim=c(n.states, n.states, Nstar, n.occasions))
for (i in 1:Nstar){
   for (t in 1:n.occasions){
      STATE_MATRIX[,,i,t] <- matrix(c(
      0,		  	0,        				0,   				       0,     #this line is not needed since it is conditional on entry state
      0,        phiA*(1-psiAB), 	phiA*psiAB,        1-phiA,
      0,        0,           			phiB, 			       1-phiB,
      0,       	0,                    		0,         1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

#Observation process
OBS_MATRIX <- array(NA, dim=c(n.states, n.obs, Nstar, n.occasions))
for (i in 1:Nstar){
   for (t in 1:n.occasions){
      OBS_MATRIX[,,i,t] <- matrix(c(
      0,  0,  1,                   #not detected if not entered
      pA, 0,  1-pA,                #detection if egg burrow
      pB*beta,  pB*(1-beta), 1-pB, #detection if chick burrow
      0,  0,  1),                  #detection if terminated
      nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Function to simulate capture-recapture data under the JS model
simul.js_ms_me <- function(STATE_MATRIX, OBS_MATRIX, Nstar, n.occasions, n.obs, n.states){

   z <- chTRUE <- matrix(0, ncol = n.occasions+1, nrow = Nstar) #observed (z) and true (true) processes including dummy occasion
   # Generate total number of entering burrows per occasion
   B <- rmultinom(n = 1, size = Nstar, prob = b) # number of entering burrows at each t
   ent.occ <- numeric()
   for (t in 1:n.occasions){
      ent.occ <- c(ent.occ, rep(t+1, B[t])) #vector of entry occasions, t at which given i entered
   }

   # Simulating survival, transition, and detection
   for (i in 1:Nstar){ #for each individual
    #prior to entry
    z[i,1:(ent.occ[i]-1)] <- 1            #all are 1s before first entry occasion
    chTRUE[i,1:(ent.occ[i]-1)] <- n.obs   #cannot be observed prior to capture

    #at entry occasion 
    for(t in ent.occ[i]){ #for t when t matches day of entry for i 
     z[i,t] <- which(rmultinom(1, 1, INIT.STATE[t-1,])==1)  #obs state at entry occasion; position of single row of init.state where i entered, which is state 2, egg
     obs <- which(rmultinom(1, 1, OBS_MATRIX[z[i,t],,i,t-1])==1) #obs given state; find the position of the 1, it is either 1 or 3, so can't enter as chick, but don't understand how/why this works; help
     chTRUE[i,t] <- obs
    }

    #after entry occasion
    if(ent.occ[i]<(n.occasions+1)){     # for occasions after entry occasion
     for(t in ent.occ[i]:n.occasions+1){
      z[i,t] <- which(rmultinom(1, 1, STATE_MATRIX[z[i,t-1],,i,t-1])==1) #state
      obs <- which(rmultinom(1, 1, OBS_MATRIX[z[i,t],,i,t-1])==1)   #obs given state
      chTRUE[i,t] <- obs
     }#t
    }#if
  }#i
   
   # Remove individuals never captured
   captured <- apply(chTRUE, 1, min) < n.obs #captured if minimum of chTRUE rows contains 1 or 2, otherwise not detected
   chOBS <- chTRUE[captured,] #keep the captured subset
   alive <- z #vector of 1s and 0s for alive/dead
    alive[z == 1|z == 4] <- 0 
    alive[z == 2|z == 3] <- 1
   Nt <- colSums(alive)    # Actual population size; add up all the 1s
   Negg <- Nchick <- NA
   for(t in 1:(n.occasions+1)){
    Negg[t] <- length(which(z[,t] == 2)) # number of eggs at any given t
    Nchick[t]<-length(which(z[,t] == 3)) # number of chicks at any given t
   }
   return(list(CH = chOBS, B = B, N = Nt, Negg = Negg, Nchick = Nchick))
} #funciton

# Execute simulation function
sim <- simul.js_ms_me(STATE_MATRIX, OBS_MATRIX, Nstar, n.occasions, n.obs, n.states)
CH <- sim$CH
 #verify
 sim$N
 sim$B
 sum(sim$B) #should equal simulated value of Nstar


# model

JS_ME_sim <- model {

# Priors and constraints
beta ~ dunif(0,1)       ## start with simple constant model
mean.psiAB ~ dunif(0,1)   ## start with simple constant model
alpha[1] ~ dunif(0,1)   ## if enter at occasion 1, probability entered as chick burrow
alphaKeep <- alpha[1]    ## this is the only one we want to keep 

for (u in 1:2){
 mean.phi[u] ~ dunif(0,1)    # Priors for mean state-spec. survival
 mean.p[u] ~ dunif(0,1)      # Priors for mean state-spec. recapture
}
for (t in 1:(n.occasions-1)){
 phiA[t]     <- mean.phi[1]   # egg survival
 phiB[t]     <- mean.phi[2]   # chick survival
 gamma[t]    ~ dunif(0,1)    # Prior for entry probabilities at occasion t
 pA[t]       <- mean.p[1]     # egg burrow detection
 pB[t]       <- mean.p[2]     # chick burrow detection
 psiAB[t]    <- mean.psiAB    # transition probability
 alpha[t+1]  <- 0              #probability of entry as chick on occasion t+1 
}

# Likelihood 
for (i in 1:M){
 # Define latent state at first occasion
 z[i,1] <- 1            # Make sure that all M individuals are in state 1 at t=1 (dummy occasion)
 for (t in 2:n.occasions){
  # State process: draw S(t) given S(t-1)
  z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
  # Observation process: draw O(t) given S(t)
  y[i,t] ~ dcat(po[z[i,t], i, t-1,])  
   #indexing example from Kery and Schaub 10.3.2, y[i,2] is a function of state[i,2] and observation array t=1
   #so if detection covariates are used, DO NOT ADD DUMMY OCCASION to detection covariates as detection data at occasion t uses array t-1

  # Derived variables
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
      po[2,i,t,1] <- pA[t]                    #probability an egg burrow is detected with a burrow visit
      po[2,i,t,2] <- 0                        #probability an egg burrow is detected with a prey visit
      po[2,i,t,3] <- 1-pA[t]                  #probability an egg burrow is not detected
      po[3,i,t,1] <- beta*pB[t]               #probability a chick burrow is detected with a burrow visit 
      po[3,i,t,2] <- (1-beta)*pB[t]           #probability a chick burrow is detected with a prey visit
      po[3,i,t,3] <- 1 - pB[t]                #probability a chick burrow is not detected
      po[4,i,t,1] <- 0                        #probability a terminated burrow is detected with a burrow visit
      po[4,i,t,2] <- 0                        #probability a terminated burrow is detected with a prey visit
      po[4,i,t,3] <- 1                        #probability a terminated burrow is not detected
  } #t
 } #M

#derive abundances (lots of options here, e.g., site specific, annual, etc.) 
 for (t in 1:(n.occasions-1)){
  N.egg[t]    <- sum(egg[1:M,t])
  N.chick[t]  <- sum(chick[1:M,t])
  N.active[t] <- sum(active[1:M,t])
  qgamma[t] <- 1-gamma[t]
  b[t] <- cprob[t] / psi      # Entry probability
}
 cprob[1] <- gamma[1]
 for (t in 2:(n.occasions-1)){
  cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
 }#t
 psi <- sum(cprob[])            # Inclusion probability

Nstar<-sum(EverActive[1:M])
 
}
 
write.model(JS_ME_sim, "JS_ME_sim.txt")
model.file = paste(getwd(),"JS_ME_sim.txt", sep="/")


# Augment data
nz <- 50
CH.aug <- rbind(CH, matrix(n.obs, nrow = nz, ncol = dim(CH)[2])) 

get.first <- function(x) min(which(x < n.obs))
f <- apply(CH, 1, get.first)
get.last <- function(x) max(which(x < n.obs))
l <- apply(CH, 1, get.last)

#inits for state process z[i,t]
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

zInitsObs <- js.me.init(CH, f, l, nz)
zInits<-rbind(zInitsObs,matrix(1,nr=nz, nc=dim(CH)[2]))
zInits[,1]<-NA #deterministic on first occasion

inits <- function(){list(mean.phi = runif(2, 0, 1), 
                         mean.p = runif(2, 0, 1), 
                         z = zInits)}    
#data
jags.data <- list(y = CH.aug, n.occasions = n.occasions+1 , M = dim(CH.aug)[1])
   
# Parameters monitored
parameters <- c("mean.phi", "mean.p", "beta", "gamma", "alphaKeep", "mean.psiAB","N.egg","N.chick", "N.active","Nstar","psi","b")

# MCMC settings
nc <- 3; nAdapt <- 500; nb <- 100;  ni <- 10000; nt <- 1

# Run
js_me <- jags(jags.data, inits, parameters, model.file = "model_JS_MS_ME.txt",
              n.chains = nc, n.adapt=nAdapt, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)

plot(js_me$samples[,"alphaKeep"])

#abundances
outmat<-as.matrix(js_me$samples)
plot(apply(outmat[,grep("N.active\\[", colnames(outmat))],2,median), pch=16, type="b", ylab="Number", xlab=c("Day of season"), ylim=c(0,Nstar))
points(apply(outmat[,grep("N.egg\\[", colnames(outmat))],2,median), pch=16, type="b", col="red")
points(apply(outmat[,grep("N.chick\\[", colnames(outmat))],2,median), pch=16, type="b", col="blue")
legend("topleft", legend=c("Active","Egg", "Chick"), text.col=c("black", "red", "blue"), bty="n")
 #true values
 points(sim$N[2:(n.occasions+1)], pch=0)
 points(sim$Negg[2:(n.occasions+1)], pch=0, col="red")
 points(sim$Nchick[2:(n.occasions+1)], pch=0, col="blue")

#estimated 'birth' rates
plot(apply(outmat[,grep("b\\[", colnames(outmat))],2,median), pch=16, type="b", ylab="Probability", xlab=c("Day of season"), ylim=c(0,1))
 points(b, pch=0) #true values




