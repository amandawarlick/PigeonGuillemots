


model_MEJS <- function () {

# Priors and constraints
for (i in 1:M) {
  for (t in 1:(n.occasions-1)){
        logit(phiA[i,t]) <- mu.phi[1] #+ b.egg.y[year_i[i]] 
        logit(phiB[i,t]) <- mu.phi[2] #+ b.chick.y[year_i[i]] 
    gamma[i,t] <- mean.gam    # Prior for entry probabilities at occasion t
    pA[i,t] <- mean.p[1]     # egg burrow detection
    pB[i,t] <- mean.p[2]     # chick burrow detection
    psiAB[i,t] <- mean.psiAB
  } #t
} #i

b ~ dunif(0,1)  # prior for assignment probability
mean.psiAB ~ dunif(0,1)   # transition from egg to chick
mean.gam ~ dunif(0,1)
pi ~ ddirch(alpha[1:3])   #dirichlet prior for multinomial
alpha[1] <- 1/3
alpha[2] <- 1/3
alpha[3] <- 1/3

for (u in 1:2){
 mu.phi[u] ~ dnorm(0, 0.001) 
 #mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
 mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
 #mean.p[u] ~ dunif(0.3, 1)      # Priors for mean state-spec. recapture

}

  b.chick.y[1] <- 0
  b.egg.y[1] <- 0

phiA.int <- 1 / (1+exp(-mu.phi[1]))
phiB.int <- 1 / (1+exp(-mu.phi[2]))

# Likelihood 
for (i in 1:M){
 # Define latent state at first occasion
  z[i,1] ~ dcat(pi[1:3])      #all M individuals have probability of being in 1 of 3 states
  #y[i,1] ~ dcat(po[z[i,1], i, 1, 1:3])  #covered below
  
 for (t in 2:n.occasions){
  # State process: draw S(t) given S(t-1); daily
  z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:4])
 } #t
   
  # Observation process: draw O(t) given S(t); n_visit approach via Nathan
  for (k in 1:visits_i[i]) {  
  y[i,k] ~ dcat(po[z[i,effort[i,k]], i, k, 1:3]) #add indices back in if modeling p
  } #k
} #i

# Define transition and observation matrices
 for (i in 1:M){
  for (t in 1:(n.occasions-1)) {
      # Define probabilities of state S(t+1) given S(t)   
      ps[1,i,t,1] <- 1-gamma[i,t]  #probability of not entering
      ps[1,i,t,2] <- gamma[i,t]                    #probability of entering as egg
      ps[1,i,t,3] <- 0                          #can't enter as a chick after day 1
      ps[1,i,t,4] <- 0                             #probability of entering as terminated
      ps[2,i,t,1] <- 0                             #probability egg goes to 'not entered'
      ps[2,i,t,2] <- (1-psiAB[i,t])*phiA[i,t]          #probability of surviving egg state and not transitioning
      ps[2,i,t,3] <- phiA[i,t]*psiAB[i,t]              #probability of surviving egg state and hatching
      ps[2,i,t,4] <- 1-phiA[i,t]                      #probability of a failed egg
      ps[3,i,t,1] <- 0                                #probability chick goes to 'not entered'
      ps[3,i,t,2] <- 0                                #probability chick goes to egg
      ps[3,i,t,3] <- phiB[i,t]
      ps[3,i,t,4] <- 1-phiB[i,t]               #probability of failed chick
      ps[4,i,t,1] <- 0                      #probability terminated goes to 'not entered'
      ps[4,i,t,2] <- 0                       #probability terminated goes to egg (maybe happens)
      ps[4,i,t,3] <- 0                       #probability terminated goes to chick
      ps[4,i,t,4] <- 1                       #probability terminated goes to terminated
  } #t

   for (t in 1:visits_i[i]) {
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0                        #'not entered' burrow is detected with a burrow visit
      po[1,i,t,2] <- 0                        #'not entered' burrow is detected with a prey visit
      po[1,i,t,3] <- 1                        #'not entered' burrow is not detected
      po[2,i,t,1] <- pA[i,t]                    #egg burrow is detected with a burrow visit
      po[2,i,t,2] <- 0                        #egg burrow is detected with a prey visit
      po[2,i,t,3] <- 1-pA[i,t]                  #egg burrow is not detected
      po[3,i,t,1] <- b * pB[i,t]                #chick burrow is detected with a burrow visit 
      po[3,i,t,2] <- (1 - b) * pB[i,t]          #chick burrow is detected with a prey visit
      po[3,i,t,3] <- 1 - pB[i,t]                #chick burrow is not detected
      po[4,i,t,1] <- 0                        #terminated burrow is detected with a burrow visit
      po[4,i,t,2] <- 0                        #terminated burrow is detected with a prey visit
      po[4,i,t,3] <- 1                        #terminated burrow is not detected
  } #t
 }#M

#mean over study period
for (i in 1:M) {
  days.chick[i] <- sum(z[i,] == 3)
  days.egg[i] <- sum(z[i,] == 2)
  fledged_high[i] <- step(days.chick[i]-32) #if step() greater than zero, 1
  fledged_low[i] <- step(days.chick[i]-45) #if step() greater than zero, 1
  everActive[i] <- max(z[i,]>1) #number of active ever
 } #i

#seasonal values
#for (i in 1:M) {
  for (t in 2:n.occasions) {
  chick[t] <- sum(z[1:M,t] == 3)
  egg[t] <- sum(z[1:M,t] == 2)
   } #t
 #} #i

#means
n.fledged.low <- sum(fledged_low[1:M])
n.fledged.high <- sum(fledged_high[1:M])
n.active.burrow <- sum(everActive[1:M])
nest.succ.low <- n.fledged.low/n.active.burrow
nest.succ.high <- n.fledged.high/n.active.burrow
mean.days.chick <- mean(days.chick)
mean.days.egg <- mean(days.egg)
} #mod

write.model(model_MEJS, "model_MEJS.txt")
model.file = paste(getwd(),"model_MEJS.txt", sep="/")


# Function to simulate capture-recapture data under a multi-event JS model;
# adapted with Nathan's help from Kery & Schaub
n.occasions <- 98            # Number of capture occasions
#Nstar <- nrow(ch_sim)     # Superpopulation size (total number of burrows active on at least one day)
Nstar <- 150
n.states <- 4                   # not entered, egg, chick, terminated
n.obs <- 3                      # burrow visit, prey delivery, not detected
gamma <- c(0.5, seq(0.2, 0.01, length = n.occasions-1))      # Entry probabilities decrease with time
gamma[(n.occasions-15):n.occasions] <- 0  #no new can enter once forced exit at t-15

phiA <- 0.99                    # survival state A (egg burrow)
phiB <- 0.98                  # survival state B (chick burrow)
psiAB <- 0.1                   # transition A to B conditional on survival
pA <- 0.5                       # detection for egg burrow
pB <- 0.8                       # detection for chick burrow
b <- 0.2                     # conditional on observing true state chick, probability obs was bv

#Define matrics

#state at entry
# init.st <- matrix(NA, nrow = n.occasions+1, ncol = n.states)
#  for(t in 1:nrow(init.st)){
#   init.st[t,]<- c(0.33, 0.33, 0.33, 0) #multinomial probability of entering into 1-3 states
#  }
init.st <- c(0, rep(1/2, 2), 0) #Equal probability of being in states 2 or 3... or 4?

#state process conidtional on entry (survival and transition)
st.proc <- array(NA, dim = c(n.states, n.states, Nstar, n.occasions))
for (i in 1:Nstar){
   for (t in 1:n.occasions){
      st.proc[,,i,t] <- matrix(c(
      (1-gamma[t]), gamma[t], 0, 0,     #this line is not needed since it is conditional on entry state?
      0, phiA*(1-psiAB), phiA*psiAB, 1-phiA,
      0, 0, phiB, 1-phiB,
      0, 0, 0, 1), nrow = n.states, byrow = TRUE)
      } #t
   } #i

#Observation process
obs.proc <- array(NA, dim=c(n.states, n.obs, Nstar, n.occasions))
for (i in 1:Nstar){
   for (t in 1:n.occasions){
      obs.proc[,,i,t] <- matrix(c(
      0,  0,  1,                   #not detected if not entered
      pA, 0,  1-pA,                #detection if egg burrow
      pB*b,  pB*(1-b), 1-pB, #detection if chick burrow
      0,  0,  1),                  #detection if terminated
      nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Function to simulate capture-recapture data under the JS model
sim.mejs <- function(st.proc, obs.proc, Nstar, n.occasions, n.obs, n.states){

   z <- ch <- matrix(0, ncol = n.occasions, nrow = Nstar)  
   B <- rmultinom(1, Nstar, gamma) #total number of entering burrows per occasion
   ent.occ <- numeric()
   for (t in 1:n.occasions){
      ent.occ <- c(ent.occ, rep(t, B[t])) #entry day/occasion for each burrow
   }

   #survival, transition, and detection
   for (i in 1:Nstar){
    #prior to entry
    z[i,1:(ent.occ[i]-1)] <- 1 #not entered until entry occasion
    ch[i,1:(ent.occ[i]-1)] <- n.obs #cannot be observed prior to capture, all 3s

    #at entry occasion 
    for (t in ent.occ[i]){ #when t matches entry day for i
     #z[i,t] <-  which(rmultinom(1, 1, init.st[t-1,])==1) #state at entry occasion
     z[i,t] <-  which(rmultinom(1, 1, init.st)==1) #state at entry occasion, either 2 or 3
     event <- which(rmultinom(1, 1, obs.proc[z[i,t],,i,t])==1) #obs given state; [states, obs, N, t)
     ch[i,t] <- event
    } #t

    #after entry occasion
    if(ent.occ[i]<(n.occasions)){ #survival process if entered prior to end of study
      for (t in ent.occ[i]:(n.occasions-1)){
      z[i,t+1] <- which(rmultinom(1, 1, st.proc[z[i,t],,i,t])==1)   #state
      event <- which(rmultinom(1, 1, obs.proc[z[i,t],,i,t])==1)     #obs given state
      ch[i,t+1] <- event
     }#t
    }#if
    
    #forced exit: should instead be forced to state 4 after X number of days either in state 2 or state 3
    z[i,c((n.occasions-15):n.occasions)] <- 4
    ch[i,c((n.occasions-15):n.occasions)] <- 3
  }#i
   
   # Remove individuals never captured
   captured <- apply(ch, 1, min) < n.obs
   chOBS <- ch[captured,]
   alive <- z
    alive[z==1|z==4] <- 0
    alive[z==2|z==3] <- 1
   Nt <- colSums(alive)    # Actual population size at time t
   Negg <- Nchick <- NA
   #for(t in 1:(n.occasions+1)){
   for(t in 1:(n.occasions)){
    Negg[t] <- length(which(z[,t]==2))
    Nchick[t] <- length(which(z[,t]==3))
   }
   return(list(CH = chOBS, B = B, N = Nt, Negg = Negg, Nchick = Nchick, z = z))
} #function

# Execute simulation function
sim <- sim.mejs(st.proc, obs.proc, Nstar, n.occasions, n.obs, n.states)
CH <- sim$CH
z <- sim$z

#verify
# sim$N
# sim$B
# sum(sim$B) #should equal simulated value of Nstar

# Augment data
nz <- (Nstar - dim(CH)[1]) + 0.2*Nstar #get back up to Nstar (capture < Nstar dep on params) + 20%
CH.aug <- rbind(CH, matrix(n.obs, nrow = nz, ncol = dim(CH)[2])) 
ch <- CH.aug
M <- dim(ch)[1]

visits_i <- rep(98/7, M) #number of visits for each nest

#survey days for time-between
effort <- matrix(NA, nrow = M, ncol = max(visits_i))
for (i in 1:M){
  #for (t in 1:n.occasions) {
    effort[i,] <- 7*(1:max(visits_i)) #even sampling 7 days apart same at all nests
  #}
}

get.first <- function(x) min(which(x < n.obs)) #first time that there is a non-3
f <- apply(CH, 1, get.first)
get.last <- function(x) max(which(x < n.obs))  #last time there is a non-3
l <- apply(CH, 1, get.last)

#inits for state process z[i,t]
js.me.init <- function(ch, f, l, nz){
  inits <- matrix(NA, nr = dim(ch)[1], nc = dim(ch)[2])
  #inits[,1] <- NA  #deterministic on occasion 1

  for(i in 1:dim(ch)[1]) { 
    inits[i,1:f[i]-1] <- 1  #1 before first capture 
    if(l[i]<dim(ch)[2]){
     inits[i, (l[i]+1):dim(ch)[2]] <- rep(3, (dim(ch)[2]-l[i])) #3 from last obs to end
    }
    if(sum(ch[i,] == 2, na.rm=T)>0){ #sites with observed prey delivery
       if(f[i]==min(which(ch[i,]  == 2))){  #if first obs is a prey delivery
        inits[i, f[i]:l[i]]<-3              #chick until end of study 
        if(f[i]>2) inits[i, 2:(f[i]-1)] <- 2  #egg prior to chick

       }else{
        inits[i, f[i]:(min(which(ch[i,]  == 2))-1)] <- 2
        inits[i, min(which(ch[i,]  == 2)):l[i]] <- 3
       }
    }else{ #never detected a prey delivery
      if(sum(ch[i,] == 2, na.rm=T)==0){      
        inits[i, f[i]:(max(which(ch[i,] == 1)))] <- 2
      }
    }
    inits[i,c((n.occasions-15):n.occasions)] <- 4 #kicked out
   }#i
  return(inits)
}

zInitsObs <- js.me.init(CH, f, l, nz)
z.init <- rbind(zInitsObs,matrix(1, nr = nz, nc = dim(CH)[2])) #add augmented

   
# Bundle data
jags.data <- list(y = ch, n.occasions = n.occasions, 
                  #nyear = length(unique(year_i)),
                  #z =z,
                  effort = effort, #from setup, survey_days
                  visits_i = visits_i, 
                  M = dim(ch)[1])

inits <- function(){list(#mean.phi = runif(2, 0, 1), 
                         #mean.p = runif(2, 0.3, 1),
                         z = z.init)}  

# Parameters monitored
parameters <- c('phiA.int', 'phiB.int', "mean.p", "b", 
                'n.active.burrow', 'n.fledged.low', 'n.fledged.high', 'mean.days.chick', 'mean.days.egg',
                "mean.psiAB", 'mu.phi',
                'nest.succ.low', 'nest.succ.high', 
                "N.active", "Nstar")
     
# MCMC settings
ni <- 3; nt <- 1; nb <- 1; nc <- 1
     
out_sim <- jags(jags.data, inits, parameters, model.file = model.file,
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)



