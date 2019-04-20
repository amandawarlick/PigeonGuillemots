
library(rjags)
library(jagsUI)
library(R2OpenBUGS)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)

counts <- read.csv('count_dat.csv', header = T, stringsAsFactors = F) 


#colony level
# count_dat <- counts %>%
#   select(year, site, week, PG_count) %>%
#   filter(week < 33 & week > 26) 
# 
# count_wide <- count_dat %>%
#   dcast(year + site ~ week, value.var = 'PG_count', fun.aggregate = mean)
# 
# cnt <- matrix(NA, nrow = dim(count_wide)[1], ncol = 8)
# for (i in 1:dim(count_wide)[1]) {
#   temp <- round(as.numeric(count_wide[i, 2 + c(which(!is.na(count_wide[i,3:dim(count_wide)[2]])))]))
#   num <- length(temp)
#   cnt[i,] <- c(temp, rep(NA, 8-num))
# }
# 
# y_cnt <- cnt[,1:3] #just start w/ 3 obs?
# site_cnt <- count_dat[,'site']
# year <- count_dat[,'year']

#island level - starting here for now
count_yr <- counts %>%
  group_by(year, week) %>%
  summarize(PG_count = sum(PG_count, na.rm = T), bv = sum(bv, na.rm = T), pv = sum(pv, na.rm = T),
            v = mean(v, na.rm = T), temp = mean(temp, na.rm = T),
            mins = mean(mins, na.rm = T), upwell = mean(upwell, na.rm = T)) 
height <- count_yr %>%
  select(year, week, v) %>%
  filter(week > 26 & week < 32) %>%
  dcast(year ~ week, value.var = 'v')
#height[11,2:6] <- NA
temp <- count_yr %>%
  select(year, week, temp) %>%
  filter(week > 26 & week < 32) %>%
  dcast(year ~ week, value.var = 'temp')
#temp[9:11,2:6] <- NA
upwell_j <- count_yr %>%
  select(year, week, upwell) %>%
  filter(week > 26 & week < 32) %>%
  dcast(year ~ week, value.var = 'upwell')
upwell_i <- count_yr %>%
  select(year, week, upwell) %>%
  group_by(year) %>%
  summarize(up = mean(upwell))

yr_cnt <- count_yr %>%
  select(year, week, PG_count) %>%
  filter(week > 26 & week < 32) %>%
  dcast(year ~ week, value.var = 'PG_count', fun.aggregate = mean)


#####null model

nMix <- function () {
  
  #priors
  lambda ~ dunif(0,10)
  p ~ dunif(0,1)

  #likelihood
  ##ecological process
  for (i in 1:R) { #R is nrow (number of years or colony-year combinations)
    N[i] ~ dpois(lambda)
    ###obs process
    for (j in 1:reps) { #reps is number of week replicates in the sample
      y[i,j] ~ dbin(p, N[i])
    } #j
    N_est[i] <- sum(N[i])
  } #i
}

write.model(nMix, "nMix.txt")
model.file = paste(getwd(),"nMix.txt", sep="/")

########run null model

y <- as.matrix(yr_cnt[,2:dim(yr_cnt)[2]])
jags.data <- list(y = y, R = nrow(y), reps = ncol(y))

Nst <- apply(y, 1, max) + 1
inits <- function(){list(N = Nst, b0 = runif(1,-1,1))}  

# Parameters monitored
parameters <- c('p', 'N_est')

# MCMC settings
ni <- 1000; nt <- 1; nb <- 500; nc <- 3

(out_null <- jags(jags.data, inits, parameters, model.file = model.file,
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T))

######covariates

nMix_covs <- function () {
  
  #priors
  for (i in 1:R) {
    log(lambda[i]) <- a0 #+ a.upwell*upwell_i[i]
    #mu.lam[i] <- exp(lambda[i]) 

    for (j in 1:reps) {
      #logit.p[i,j] <- b0 #+ b.upwell*upwell_j[i,j]
      logit(logit.p[i,j]) <- b0 #+ b.upwell*upwell_j[i,j]  #plogit on logit scale
      p[i,j] <- exp(logit.p[i,j])/(1+exp(logit.p[i,j])) #equivalent to plogis(); p must be on probability scale
    }
  }
  
  b0 ~ dunif(-10, 10)
  #b.upwell ~ dunif(-10,10)
  a0 ~ dunif(-10,10)
  #a.upwell ~ dunif(-10,10) #effect on N
  
  #likelihood
  ##ecological process
  for (i in 1:R) { #R is nrow (number of years or colony-year combinations)
    N[i] ~ dpois(lambda[i])
    ###obs process
    for (j in 1:reps) { #reps is number of week replicates in the sample
      y[i,j] ~ dbin(p[i,j], N[i])
    } #j
    N_est[i] <- sum(N[i])
  } #i
}

write.model(nMix_covs, "nMix_covs.txt")
model.file = paste(getwd(),"nMix_covs.txt", sep="/")

#### run covariate model

# Bundle data
y <- as.matrix(yr_cnt[,2:dim(yr_cnt)[2]])
jags.data <- list(y = y, 
                  temp = as.matrix(temp[,2:6]), 
                  upwell_j = as.matrix(upwell_j[,2:6]), 
                  upwell_i = as.matrix(upwell_i[,2]),
                  height = as.matrix(height[,2:6]),
                  R = nrow(y), reps = ncol(y))

Nst <- apply(y, 1, max) + 1
inits <- function(){list(N = Nst, b0 = runif(1,-1,1))}  

# Parameters monitored
parameters <- c('p.mu', 'p.mean', 'b0', 'a0', 'N_est', 'b.upwell', 'b.height', 'b.up.p', 'lambda')

# MCMC settings
ni <- 1000; nt <- 1; nb <- 500; nc <- 3

(out <- jags(jags.data, inits, parameters, model.file = model.file,
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T))


