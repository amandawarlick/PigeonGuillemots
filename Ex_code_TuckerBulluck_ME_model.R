# JAGS model for estimating annual survival of Prothonotary Warblers 
# and testing for an effect of conspecific brood parasitism on host
# females

# Anna Tucker (annamtucker@gmail.com)

# Tucker, A.M. and L.P. Bulluck. 2017. No evidence for a negative 
# effect of conspecific brood parasitism on Prothonotary Warbler
# females. Ibis.


model {
  
  # for shrinkage of uninformative priors on beta
  betasig ~ dgamma(1, 0.1)
  betatau <- 1/betasig^2
  
  # for indicator variables and beta coefficients
  # 1 = host
  # 2 = age
  # 3 = nestlings
  
  for(i in 1:3){
    ind[i] ~ dbern(0.5)
    b[i] ~ dnorm(0, betatau)
    beta[i] <- ind[i]*b[i]
  }
  
  # set up state and observation matrices
  for(i in 1:N){
    for(t in f[i]:(n.years-1)){
      
      logit(phiN[i,t]) <- mu.s + beta[2]*age[i,t] + 
        beta[3]*nestlings[i,t] + 
        epsilon[t]
      logit(phiH[i,t]) <- mu.s + beta[1] + 
        beta[2]*age[i,t] + 
        beta[3]*nestlings[i,t] + 
        epsilon[t]
    }
  }
  
  for(t in 1:(n.years-1)){
    epsilon[t] ~ dnorm(0, tau) #annual variation error term
  }
  
  mu.s <- log(mean.phi/(1-mean.phi))
  mean.phi ~ dunif(0, 1) #baseline survival
  
  logit(mean.phiN) <- mu.s
  logit(mean.phiH[1]) <- mu.s + beta[1]
  
  # initial state assignment
  init[1] <- pi
  init[2] <- 1-pi
  init[3] <- 0
  
  # observation process
  
  OBS[1,1] <- 1-p
  OBS[1,2] <- p*delta
  OBS[1,3] <- 0
  OBS[1,4] <- p*(1-delta)
  
  OBS[2,1] <- 1-p
  OBS[2,2] <- 0
  OBS[2,3] <- p*delta
  OBS[2,4] <- p*(1-delta)
  
  OBS[3,1] <- 1
  OBS[3,2] <- 0
  OBS[3,3] <- 0
  OBS[3,4] <- 0
  
  OBS.init[1,1] <- 0
  OBS.init[1,2] <- delta
  OBS.init[1,3] <- 0
  OBS.init[1,4] <- 1-delta
  
  OBS.init[2,1] <- 0
  OBS.init[2,2] <- 0
  OBS.init[2,3] <- delta
  OBS.init[2,4] <- 1-delta
  
  OBS.init[3,1] <- 1
  OBS.init[3,2] <- 0
  OBS.init[3,3] <- 0
  OBS.init[3,4] <- 0
  
  # state process
  for(i in 1:N){
    for(t in f[i]:(n.years-1)){
      
      STATE[1,i,t,1] <- phiN[i,t]*(1-psiNH)
      STATE[1,i,t,2] <- phiN[i,t]*psiNH
      STATE[1,i,t,3] <- 1-phiN[i,t]
      
      STATE[2,i,t,1] <- phiH[i,t]*psiHN
      STATE[2,i,t,2] <- phiH[i,t]*(1-psiHN)
      STATE[2,i,t,3] <- 1-phiH[i,t]
      
      STATE[3,i,t,1] <- 0
      STATE[3,i,t,2] <- 0
      STATE[3,i,t,3] <- 1
    }
  }
  
  # likelihood
  for (i in 1:N){
    z[i,f[i]] ~ dcat(init[1:3])
    y[i,f[i]] ~ dcat(OBS.init[z[i,f[i]],1:4])
    
    for (j in (f[i]+1):n.years){
      z[i,j] ~ dcat(STATE[z[i,j-1],i,j-1,1:3])
      y[i,j] ~ dcat(OBS[z[i,j],1:4])
    }
  }
  
  # prior
  sigma <- sqrt(1/tau)
  tau ~ dgamma(0.1, 0.1)
  
  pi ~ dunif(0, 1)
  p ~ dunif(0,1)
  delta ~ dunif(0,1)
  psiNH ~ dunif(0, 1)
  psiHN ~ dunif(0, 1)
}

params = c("p", "b", "pi", "delta", "psiNH", "psiHN",
           "mean.phi", "sigma", "mean.phiH", "betasig",
           "ind", "mean.phiN")


## generate initial values for latent states
# z1 - assumes all unknowns are nonhosts
# z2 - assumes all unknowns are hosts
z = CH
z [z==3] = 2
for (i in 1:N) {
  for (j in 1:years) {
    if (j > f[i] & CH[i,j]==0) {z[i,j] = which(rmultinom(1, 1, c(0.5, 0.5))==1)}
    if (j < f[i]) {z[i,j] = NA}
  }
}
z1 = as.matrix(z)

z = CH
z [z==3] = 1
for (i in 1:N) {
  for (j in 1:years) {
    if (j > f[i] & CH[i,j]==0) {z[i,j] = which(rmultinom(1, 1, c(0.5, 0.5))==1)}
    if (j < f[i]) {z[i,j] = NA}
  }
}
z2 = as.matrix(z)

init1 = list(z=z1)
init2 = list(z=z2)
inits = list(init1, init2, init1, init2)


