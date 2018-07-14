################################################
#
# First steps towards integrated population models
#
# IPM workshop, Aberdeen June 2018
#
################################################

# Michael Schaub, 1 June 2018


# Load library
library(jagsUI)


# Specify path
setwd('...')

# Load shrike data

data <- load("WoodchatShrike.Rdata")



#############################################


# 5.2. Use of demographic data in the analysis of a population projection model

# 5.2.1. Combining capture-recapture data with a population projection model

# Specify the model in BUGS language
cat(file = "m1.jags", "
    model { 
    # Priors and constraints
    mean.sj ~ dunif(0, 1)
    mean.sa ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)
    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]            # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    
    # Population model
    # Model for initial state
    N[1,1] <- 1
    N[2,1] <- 1
    
    # Loop over time
    for (t in 1:T){
    # Population projection
    N[1,t+1] <- f * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
    # Annual growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])    
    }
    lambda <- ann.growth.rate[T]
    }
    ")


# Bundle data
bugs.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = dim(marray)[2], rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]), f = 1.6, T = 20)

# Initial values
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "lambda")


# MCMC settings
ni <- 3000; nt <- 1; nb <- 1000; nc <- 3

# Call JAGS from R (jagsUI)
m1 <- jags(bugs.data, inits, parameters, "m1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(m1, 3)


#############################################

# 5.2.2. Combining capture-recapture and productivity data with a population projection model

# Specify the model in BUGS language
cat(file = "m2.jags", "
    model { 
    # Priors and constraints
    mean.sj ~ dunif(0, 1)
    mean.sa ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)
    mean.f ~ dunif(0, 10)
    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    
    # Poisson regression model for productivity data
    for (i in 1:n.J){
    J[i] ~ dpois(mean.f)
    }
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]            # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    
    # Population model
    # Model for initial state
    N[1,1] <- 1
    N[2,1] <- 1
    
    # Loop over time
    for (t in 1:T){
    # Population projection
    N[1,t+1] <- mean.f * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
    # Annual growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])    
    }
    lambda <- ann.growth.rate[T]
    }
    ")

# Bundle data
bugs.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = dim(marray)[2], rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]), J = J, n.J = length(J), T = 20)

# Initial values
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p",  "mean.f", "lambda")

# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS from R (jagsUI)
m2 <- jags(bugs.data, inits, parameters, "m2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(m2, 3)


######################################

# 5.3. The first integrated population model

# Specify the model in BUGS language
cat(file = "m3.jags", "
    model { 
    # Priors and constraints
    mean.sj ~ dunif(0, 1)
    mean.sa ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)
    mean.f ~ dunif(0, 10)
    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    
    sigma.obs ~ dunif(0.5, 50)
    tau.obs <- pow(sigma.obs, -2)
    
    # State-space model for count data
    # Model for the initial population size: discrete uniform priors
    N[1,1] ~ dunif(1, 300)
    N[2,1] ~ dunif(1, 300)
    
    # Process model over time
    for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
    }
    
    # Observation model
    for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
    }
    
    # Poisson regression model for productivity data
    for (i in 1:n.J){
    J[i] ~ dpois(mean.f)
    }
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    
    # Derived parameters
    # Annual population growth rate
    for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])    
    }
    # Total population size
    for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    }
    }
    ")


# Bundle data
bugs.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = dim(marray)[2], rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]), J = J, n.J = length(J), count = count)


# Initial values
inits <- function(){list(mean.sj = runif(1, 0, 0.5), mean.sa = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 12000; nt <- 6; nb <- 2000; nc <- 3

# Call JAGS from R (jagsUI)
m3 <- jags(bugs.data, inits, parameters, "m3.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(m3, 3)



u <- col2rgb("grey92")
col.pol <- rgb(u[1], u[2], u[3], alpha = 75, maxColorValue = 255)
par(cex = 1.2)
plot(m3$mean$Ntot, type = "n", ylim = range(c(m3$q2.5$Ntot, m3$q97.5$Ntot)), ylab = "Population size", xlab = "Year", las = 1, cex = 1.5)
points(m3$mean$Ntot, type = "b", col = "black", pch = 16, lty = 1, cex = 1.5)
points(count, type = "b", col = "cornflowerblue", pch = 1, lty = 2, cex = 1.5)
points(ind$Nu["Total",], type = "b", col = "red", pch = 1, lty = 2, cex = 1.5)
T <- length(m3$mean$Ntot)
polygon(c(1:T, T:1), c(m3$q2.5$Ntot, m3$q97.5$Ntot[T:1]), border = NA, col = col.pol)
legend("topleft", legend = c("Truth", "Counts", "Estimates"), pch = c(1, 1, 16), col = c("red", "cornflowerblue", "black"), lty = c(2, 2, 1), bty = "n")




##################################

# Exercises

##################################

# Exercise 1: 
# Fit an IPM that includes demographic stochasticity


###################

# Exercise 2: 
# Fit an IPM in which adult survival is variable over time and juvenile survival shows a linear temporal trend


############################

# Exercise 3: 
# Fit an IPM that uses the Poisson distribution for the observation model of the state-space model



#######################

# Exercise 4: 
# Fit an IPM to the shrike data. Assume that only capture-recapture data of the adults [CH.A.marray], counts and productivity data are available.


######################


# Exercise 5: 
# Fit an IPM to the woodchat shrike data. Assume now, that the shikes start to breed only when they are 2 years old. 



##########################

# Exercise 6: 
# Fit an IPM to the woodchat shrike data. Still use the pre-breeding census model, but change the model in such a way that we monitor in addition the number of fledglings.











########################################

# 5.5. Estimation of demographic parameters without explicit data

# Specify the model in BUGS language
cat(file = "m4.jags", "
    model { 
    # Priors and constraints
    mean.sj ~ dunif(0, 1)
    mean.sa ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)
    mean.f ~ dunif(0, 10)
    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    
    sigma.obs ~ dunif(0.5, 50)
    tau.obs <- pow(sigma.obs, -2)
    
    # State-space model for count data
    # Model for the initial population size
    N[1,1] ~ dunif(1, 300)
    N[2,1] ~ dunif(1, 300)
    
    # Process model over time
    for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
    }
    
    # Observation model
    for (t in 1:n.occasions){
    count[t] ~ dnorm(N[1,t] + N[2,t], tau.obs)
    }
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    
    # Derived parameters
    # Annual population growth rate
    for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])    
    }
    # Total population size
    for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    }
    }
    ")

# Bundle data
bugs.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = dim(marray)[2], rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]), count = count)

# Initial values
inits <- function(){list(mean.sj = runif(1, 0, 0.5), mean.sa = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 12000; nt <- 6; nb <- 2000; nc <- 3

# Call JAGS from R (jagsUI)
m4 <- jags(bugs.data, inits, parameters, "m4.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(m4, 3)


