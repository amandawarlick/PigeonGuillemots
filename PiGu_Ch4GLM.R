library(tidyr)
library(data.table)
library(ggfortify)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
# library(stringr)
# library(magrittr) 
#library(ggTimeSeries)
library(stats) 
library(zoo)
library(sciplot)
library(R2OpenBUGS)
WINE="/usr/local/Cellar/wine/3.0_1/bin/wine"
WINEPATH="/usr/local/Cellar/wine/3.0_1/bin/winepath"
OpenBUGS.pgm="/Users/amandawarlick/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
setwd("~/Documents/R/SAFS/PigeonGuillemots")

#Ch. 3: Basicl Poisson & Binomial GLM

####### Simulated data ##############
#linear predictor is cubic binomial function
#data to go back to if trying to figure out if something works

# data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014) {
#   #generating values of time covariate
#   year <- 1:n 
#   #signal part of GLM
#   log.expected.count <- alpha + beta1*year + beta2*year^2 + beta3*year^3 #Q:why is this additive and not just yr^3?
#   expected.count <- exp(log.expected.count)
#   #poisson noise around expected counts
#   C <- rpois(n = n, lambda = expected.count) 
#   #plot simulated data; Q: dots are supposed to be "observed data", which is the poisson noise C
#   plot(year, C, type = "b", col = "black", las = 1)
#   lines(year, expected.count, col = "red")
#   return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3,
#               year = year, expected.count = expected.count, C = C))
# }
# data <- data.fn()

# mean.year <- mean(1:length(data$year))
# sd.year <- sd(1:length(data$year))
# win.data <- list(C = data$C, n = length(data$C), year = (1:length(data$year) - mean.year)/sd.year)

#analyze sim data using R

# fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
# summary(fm)


##################### Use real data ##################

#Basic Poisson GLM for (a) adult counts, (b) fecundity (successful fledglings), and 
# (c) binomial GLM for proportion of successful nesting pairs

#PiGu data
PG_data_BUGS <- read.csv("PG_data_all.csv", header = T) %>%
  select(year, date, site, week_study, PG_count, intern_data) %>%
  filter(intern_data != 'Y') %>%
  distinct() %>%
  filter(site %in% c("Cliffside", "Double Bluff North", "Double Bluff South", "Forbes Point", "Fort Casey",
                     "Harrington North", "Harrington South", "Hastie Lake South", "Keystone", "Lagoon North #1",
                     "Lagoon North #2", "Lagoon North #3", "Lagoon South", "Ledgewood", "Malmo Bluff", 
                     "Maylor Point", "Monroe Landing", "Mutiny Sands", "Possession Point", "Pratts Bluff",
                     "Rolling Hills #1", "Rolling Hills #2", "Shore Meadows", "Swantown"))

PG_yearly_island <- PG_data_BUGS %>%
  filter(!is.na(PG_count)) %>%
  group_by(site, year) %>%
  summarize(mean_cnt = mean(PG_count)) %>%
  group_by(year) %>%
  summarize(PG_count = sum(mean_cnt))

data <- PG_yearly_island

### (a): poisson GLM for adult counts

#bundle up the data
attach(data)
mean.year <- mean(1:length(year))
sd.year <- sd(1:length(year))
win.data <- list(C = PG_count, n = length(PG_count),
                 year = (1:length(year) - mean.year)/sd.year) #center/standardize year covariate

#specify model in BUGS
GLM_Poisson_PG <- function() {
  
  #priors
  alpha ~ dunif(-20, 20)
  beta1 ~ dunif(-10, 10)
  beta2 ~ dunif(-10, 10)
  beta3 ~ dunif(-10, 10)
  #likelihood
  for(i in 1:n) {
    C[i] ~ dpois(lambda[i]) #distribution of random element
    log(lambda[i]) <- log.lambda[i] #link function
    log.lambda[i] <- alpha + beta1*year[i] + 
      beta2*pow(year[i], 2) + beta3*pow(year[i], 3) #linear predictor
  } #i
}

write.model(GLM_Poisson_PG, "GLM_Poisson_PG.txt")
model.filePG = paste(getwd(),"GLM_Poisson_PG.txt", sep="/")

#initial values - must define 1
inits <- function () list(alpha = runif(1, -2, 2),
                          beta1 = runif(1, -3, 3)) 

#parameters
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

#MCMC
ni <- 1000
nt <- 2
nb <- 50
nc <- 3

out <- bugs(data = win.data, inits = inits, 
            parameters.to.save = params,
            model.file = model.filePG, n.chains = nc, 
            n.thin = nt, n.iter = ni,
            n.burnin = nb, OpenBUGS.pgm=OpenBUGS.pgm,
            WINE=WINE,
            WINEPATH=WINEPATH,
            useWINE=TRUE,
            debug = TRUE)

print(out, dig = 3)
# mean     sd    2.5%     25%     50%     75%   97.5%  Rhat n.eff
# alpha        6.323  0.021   6.285   6.309   6.322   6.336   6.364 1.001  2800
# beta1       -0.130  0.037  -0.201  -0.156  -0.131  -0.105  -0.056 1.008   270
# beta2       -0.043  0.017  -0.077  -0.054  -0.043  -0.031  -0.009 1.001  2800
# beta3        0.014  0.022  -0.029  -0.001   0.014   0.029   0.055 1.012   180
# lambda[1]  588.029 21.973 546.845 572.600 587.600 602.800 631.877 1.006   360
# lambda[2]  598.619 13.568 571.700 589.600 598.800 607.700 624.700 1.001  2800
#...
# lambda[10] 437.319 18.594 402.900 424.300 436.700 449.800 475.600 1.005   460
# deviance   101.497  2.870  97.990  99.442 100.800 102.900 108.777 1.001  2800
# 
# DIC info (using the rule, pD = Dbar-Dhat)
# pD = 3.971 and DIC = 105.500

wingBUGS.predictions <- out$mean$lambda
plot(year, PG_count, type = "b", lwd = 1, xlab = "", ylab = "Adult Count", pch = 20)
lines(year, wingBUGS.predictions, type = "l", lwd = 1, lty = 2)

## (b) modeling fecundity

## (c) modeling binomial nesting success


############## ~~~~ ############# ~~~~ ###############
############## ~~~~ ############# ~~~~ ###############
############## ~~~~ ############# ~~~~ ###############

#Ch. 4: GLMMS - adding fixed and random effects

#(a) random year effect

PG_weekly_col <- PG_data_BUGS %>%
  select(-c(date, intern_data))

#data <- PG_weekly_col

attach(data)
mean.year <- mean(1:length(year))
sd.year <- sd(1:length(year))
win.data <- list(C = PG_count, n = length(PG_count),
                 year = (1:length(year) - mean.year)/sd.year) #center/standardize year covariate

#specify model in BUGS
GLMM_Poisson_PG <- function() {
  
  #priors
  alpha ~ dunif(-20, 20)
  beta1 ~ dunif(-10, 10)
  beta2 ~ dunif(-10, 10)
  beta3 ~ dunif(-10, 10)
  tau <- 1/(sd*sd)
  sd ~ dunif(0, 3)
  
  #likelihood
  for(i in 1:n) {
    C[i] ~ dpois(lambda[i]) #distribution of random element
    log(lambda[i]) <- log.lambda[i] #link function
    log.lambda[i] <- alpha + beta1*year[i] + beta2*pow(year[i], 2) + 
      beta3*pow(year[i], 3) + eps[i] #linear predictor & random year effect
    eps[i] ~ dnorm(0, tau) #defines random effects distribution
  } 
}

write.model(GLMM_Poisson_PG, "GLMM_Poisson_PG.txt")
model.filePG = paste(getwd(),"GLMM_Poisson_PG.txt", sep="/")

#initial values - m ust define 1
inits <- function () list(alpha = runif(1, -2, 2),
                          beta1 = runif(1, -3, 3),
                          sd = runif(1, 0, 1)) 

#parameters
params <- c("alpha", "beta1", "beta2", "beta3", "lambda", "sd", "eps")

#MCMC
ni <- 100000
nt <- 2
nb <- 5000
nc <- 3

out <- bugs(data = win.data, inits = inits, 
            parameters.to.save = params,
            model.file = model.filePG, n.chains = nc, 
            n.thin = nt, n.iter = ni,
            n.burnin = nb, useWINE=TRUE, OpenBUGS.pgm=OpenBUGS.pgm,
            WINE=WINE, WINEPATH=WINEPATH,
            debug = TRUE)

print(out, dig = 3)
# mean     sd    2.5%     25%     50%     75%   97.5%  Rhat  n.eff
# alpha        3.401  1.382   0.040   2.683   3.487   4.338   5.699 1.275     13
# beta1       -1.418  2.965  -7.513  -3.459  -1.268   0.646   4.014 1.693      6
# beta2       -1.664  1.094  -3.587  -2.416  -1.774  -0.929   0.517 1.076     33
# beta3        0.641  1.648  -2.195  -0.553   0.549   1.742   4.130 1.615      7
# lambda[1]  617.925 24.851 570.200 601.000 617.600 634.500 667.500 1.001  97000
# lambda[2]  549.257 23.422 504.200 533.400 548.900 564.800 596.200 1.001  84000
# lambda[3]  588.740 24.287 542.200 572.200 588.400 604.900 637.500 1.001 280000
#...
# lambda[8]  433.232 20.784 393.400 419.000 432.900 447.100 474.800 1.001  40000
# lambda[9]  464.248 21.486 423.000 449.600 463.900 478.600 507.400 1.001  59000
# lambda[10] 451.886 21.245 411.100 437.400 451.600 466.000 494.400 1.001 120000
# sd           2.841  0.180   2.345   2.790   2.899   2.959   2.996 1.069    120
# eps[1]       6.698  2.790   1.898   4.586   6.569   8.507  12.430 1.011    230
# eps[2]       4.482  1.779   0.957   3.267   4.380   5.697   7.992 1.366      9
#...
# eps[9]       5.612  1.697   2.175   4.520   5.596   6.742   8.899 1.071     33
# eps[10]      6.391  2.416   1.876   4.659   6.431   8.123  10.910 1.128     21
# deviance    91.211  4.467  84.460  87.960  90.560  93.780 101.700 1.001 130000

# DIC info (using the rule, pD = Dbar-Dhat)
# pD = 9.993 and DIC = 101.200


#(b) variability among groups - year and site random effects

#look at tit.txt to see layouts
# tits <- read.table("tits.txt", header = T)
# C_tits <- as.matrix(tits[5:13])

PG_weekly_col <- PG_data_BUGS %>%
  select(-c(date, intern_data))

PG_col_wide <- PG_weekly_col %>%
  group_by(site, year) %>%
  summarize(PG_count = max(PG_count)) %>%
  dcast(site ~ year, value.var = 'PG_count')

#adjust when add new years of data
C <- as.matrix(PG_col_wide[2:11])

#start with intercept-only model - expect constant over space/time
GLMM_0 <- function() {
  #priors
  alpha ~ dnorm(0, 0.01) #log (mean count)
  #likelihood
  for (i in 1:nyear) {
    for (j in 1:nsite) {
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha
    } #j
  } #i
}
write.model(GLMM_0, "GLMM0.txt")
model.fileGLMM0 = paste(getwd(),"GLMM0.txt", sep="/")

#bundle etc.
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

inits <- function() list(alpha = runif(1, -10, 10))
params <- c("alpha")

ni <- 10000
nt <- 2
nb <- 100
nc <- 3

out_0 <- bugs(data = win.data, inits = inits, 
              parameters.to.save = params,
              model.file = model.fileGLMM0, n.chains = nc, 
              n.thin = nt, n.iter = ni,
              n.burnin = nb, useWINE=TRUE, OpenBUGS.pgm=OpenBUGS.pgm,
              WINE=WINE, WINEPATH=WINEPATH,
              debug = TRUE) 

print(out_0, 3)
# Current: 3 chains, each with 10000 iterations (first 100 discarded), n.thin = 2
# Cumulative: n.sims = 29700 iterations saved
# mean    sd     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# alpha       3.715 0.012    3.693    3.707    3.716    3.723    3.739 1.003  1000
# deviance 3078.771 1.405 3078.000 3078.000 3078.000 3079.000 3083.000 1.007  1100
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = Dbar-Dhat)
# pD = 1.017 and DIC = 3080.000

###############
#fixed site effects (one-way ANOVA for poisson response)
GLMM_site_f <- function() {
  #priors
  for (j in 1:nsite) {
    alpha[j] ~ dnorm(0, 0.01) #site effect
  }
  #likelihood
  for (i in 1:nyear) {
    for (j in 1:nsite) {
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
    } #j
  } #i
}
write.model(GLMM_site_f, "GLMM_site_f.txt")
model.fileGLMM_site_f = paste(getwd(),"GLMM_site_f.txt", sep="/")

#bundle etc.
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

inits <- function() list(alpha = runif(235, -1, 1))
params <- c("alpha")

ni <- 10000
nt <- 2
nb <- 1000
nc <- 3

#not working
out_site_f <- bugs(win.data, inits, params, model.file = model.fileGLMM_site_f, 
                   n.chains = nc, n.thin = nt, n.iter = ni,
                   n.burnin = nb, debug = TRUE, OpenBUGS.pgm = OpenBUGS.pgm, 
                   WINE = WINE, WINEPATH = WINEPATH, useWINE = TRUE) 

print(out_site_f, 3)

#############
#fixed site and fixed year (two-way main effects ANOVA for Poisson response)
GLMM_siteyr_f <- function() {
  #priors
  for (j in 1:nsite) {
    alpha[j] ~ dnorm(0, 0.01) #site effect
  }
  for (i in 2:nyear) {
    eps[i] ~ dnorm(0, 0.01)
  }
  eps[1] <- 0 #aliased; set value of first level of year effects to zero to avoid overparameterization?
  #likelihood
  for (i in 1:nyear) {
    for (j in 1:nsite) {
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j] + eps[i]
    } #j
  } #i
}
write.model(GLMM_siteyr_f, "GLMM_siteyr_f.txt")
model.fileGLMM_siteyr_f = paste(getwd(),"GLMM_siteyr_f.txt", sep="/")

#bundle etc.
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

inits <- function() list(alpha = runif(235, -1, 1), 
                         eps = c(NA, runif(8, -1, 1)))
params <- c("alpha", "eps")

ni <- 12000
nt <- 2
nb <- 10
nc <- 3

#not working
out_siteyr_f <- bugs(win.data, inits, params, n.chains = nc, n.thin = nt, n.iter = ni,
                     n.burnin = nb, debug = T, OpenBUGS.pgm = OpenBUGS.pgm, 
                     model.file = model.fileGLMM_siteyr_f, 
                     WINE = WINE, WINEPATH = WINEPATH, useWINE = TRUE) 

print(out_siteyr_f, 3) 

################
#random effects - site only

GLMM_site_r <- function() {
  #priors
  for (j in 1:nsite) {
    alpha[j] ~ dnorm(mu.alpha, tau.alpha) #random site effects
  }
  mu.alpha ~ dnorm(0, 0.01)
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 5)
  #likelihood
  for (i in 1:nyear) {
    for (j in 1:nsite) {
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
    } #j
  } #i
}
write.model(GLMM_site_r, "GLMM_site_r.txt")
model.fileGLMM_site_r = paste(getwd(),"GLMM_site_r.txt", sep="/")

#bundle etc.
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

inits <- function() list(mu.alpha = runif(1, 2, 3))
params <- c("alpha", "mu.alpha", "sd.alpha")

ni <- 10000
nt <- 2
nb <- 100
nc <- 3

out_site_r <- bugs(data = win.data, inits = inits, 
                   parameters.to.save = params,
                   model.file = model.fileGLMM_site_r, n.chains = nc, 
                   n.thin = nt, n.iter = ni,
                   n.burnin = nb, useWINE=TRUE, OpenBUGS.pgm=OpenBUGS.pgm,
                   WINE=WINE, WINEPATH=WINEPATH,
                   debug = TRUE) 

print(out_site_r, 3) #converges
# Current: 3 chains, each with 10000 iterations (first 100 discarded), n.thin = 2
# Cumulative: n.sims = 29700 iterations saved
# mean    sd     2.5%      25%      50%      75%    97.5%  Rhat n.eff
# alpha[1]     3.592 0.054    3.483    3.556    3.593    3.630    3.696 1.001  9000
# alpha[2]     4.445 0.037    4.373    4.421    4.446    4.470    4.516 1.001  4500
# alpha[3]     2.840 0.085    2.671    2.783    2.842    2.899    3.000 1.002  3400
# ...
# mu.alpha     3.593 0.101    3.394    3.527    3.593    3.659    3.795 1.001 30000
# sd.alpha     0.485 0.079    0.359    0.429    0.476    0.530    0.666 1.001 30000
# deviance  1566.402 7.014 1555.000 1561.000 1566.000 1571.000 1582.000 1.001 29000
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = Dbar-Dhat)
# pD = 23.650 and DIC = 1590.000


################
#random effects - site and year

GLMM_siteyr_r <- function() {
  #priors
  mu ~ dnorm(0, 0.01) #grand mean
  for (j in 1:nsite) {
    alpha[j] ~ dnorm(0, tau.alpha) #random site effects
  }
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 5)
  
  for (i in 1:nyear) {
    eps[i] ~ dnorm(0, tau.eps) #random year effects
  }
  tau.eps <- 1/(sd.eps*sd.eps)
  sd.eps ~ dunif(0,3)
  
  #likelihood
  for (i in 1:nyear) {
    for (j in 1:nsite) {
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + alpha[j] + eps[i]
    } #j
  } #i
}
write.model(GLMM_siteyr_r, "GLMM_siteyr_r.txt")
model.fileGLMM_siteyr_r = paste(getwd(),"GLMM_siteyr_r.txt", sep="/")

#bundle etc.
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

inits <- function() list(mu = runif(1, 0, 4), 
                         alpha = runif(235, -1, 1),
                         eps = runif(9, -1, 1))
params <- c("mu", "alpha", "eps", "sd.alpha", "sd.eps")

ni <- 10000
nt <- 2
nb <- 100
nc <- 3

out_siteyr_r <- bugs(data = win.data, inits = inits, 
                     parameters.to.save = params,
                     model.file = model.fileGLMM_siteyr_r, n.chains = nc, 
                     n.thin = nt, n.iter = ni,
                     n.burnin = nb, useWINE=TRUE, OpenBUGS.pgm=OpenBUGS.pgm,
                     WINE=WINE, WINEPATH=WINEPATH,
                     debug = TRUE) 

#doesn't converge
print(out_siteyr_r, 3)

####################

# random site/year effects and overall linear time trend

GLMM_siteyr_r_trend <- function() {
  #priors
  mu ~ dnorm(0, 0.01) #overall intercept
  beta1 ~ dnorm(0, 0.01) #overall trend
  
  for (j in 1:nsite) {
    alpha[j] ~ dnorm(0, tau.alpha) #random site effects
  }
  tau.alpha <- 1/(sd.alpha*sd.alpha)
  sd.alpha ~ dunif(0, 3)
  
  for (i in 1:nyear) {
    eps[i] ~ dnorm(0, tau.eps) #random year effects
  }
  tau.eps <- 1/(sd.eps*sd.eps)
  sd.eps ~ dunif(0, 1)
  
  #likelihood
  for (i in 1:nyear) {
    for (j in 1:nsite) {
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1*year[i] + alpha[j] + eps[i]
    } #j
  } #i
}
write.model(GLMM_siteyr_r_trend, "GLMM_siteyr_r_trend.txt")
model.fileGLMM_siteyr_r_trend = paste(getwd(),"GLMM_siteyr_r_trend.txt", sep="/")

#bundle etc.
mean.year <- mean(1:length(PG_yearly_island$year))
sd.year <- sd(1:length(PG_yearly_island$year))

win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), 
                 year = (1:length(PG_yearly_island$year) - mean.year)/sd.year)

inits <- function() list(mu = runif(1, 0, 4), 
                         alpha = runif(235, -1, 1),
                         beta1 = runif(1, -1, 1),
                         eps = runif(9, -1, 1))
params <- c("mu", "beta1", "alpha", "eps", "sd.alpha", "sd.eps")

ni <- 12000
nt <- 6
nb <- 6000
nc <- 3

out_siteyr_r_trend <- bugs(data = win.data, inits = inits, 
                     parameters.to.save = params,
                     model.file = model.fileGLMM_siteyr_r_trend, n.chains = nc, 
                     n.thin = nt, n.iter = ni,
                     n.burnin = nb, useWINE=TRUE, OpenBUGS.pgm=OpenBUGS.pgm,
                     WINE=WINE, WINEPATH=WINEPATH,
                     debug = TRUE) 
#doesn't converge
print(out_siteyr_r_trend, 3)



