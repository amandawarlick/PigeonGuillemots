
library(R2OpenBUGS)
library(dplyr)
WINE="/usr/local/Cellar/wine/3.0_1/bin/wine"
WINEPATH="/usr/local/Cellar/wine/3.0_1/bin/winepath"
OpenBUGS.pgm="/Users/amandawarlick/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

PG_weekly_colony <- PG_data_BUGS %>%
  select(-c(date, intern_data))

PG_col_wide <- PG_weekly_colony %>%
  group_by(site, year) %>%
  summarize(PG_count = max(PG_count)) %>%
  dcast(site ~ year, value.var = 'PG_count')

C <- as.matrix(PG_col_wide[2:10])

# (a) start with intercept-only model - expect constant over space/time
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
# DIC is an estimate of expected predictive error (lower deviance is better).

###############
#(b) fixed site effects (one-way ANOVA for poisson response)
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

ni <- 1000
nt <- 2
nb <- 100
nc <- 3

#not working
out_site_f <- bugs(win.data, inits, params, model.file = model.fileGLMM_site_f, 
                   n.chains = nc, n.thin = nt, n.iter = ni,
                   n.burnin = nb, debug = TRUE, OpenBUGS.pgm = OpenBUGS.pgm, 
                   WINE = WINE, WINEPATH = WINEPATH, useWINE = TRUE) 

print(out_site_f, 3)

#############
#(c) fixed site and fixed year (two-way main effects ANOVA for Poisson response)
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
#(d) random effects - site only

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

print(out_site_r, 3)
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
# DIC is an estimate of expected predictive error (lower deviance is better).


################
#(e) random effects - site and year

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

#not working
out_siteyr_r <- bugs(data = win.data, inits = inits, 
                     parameters.to.save = params,
                     model.file = model.fileGLMM_siteyr_r, n.chains = nc, 
                     n.thin = nt, n.iter = ni,
                     n.burnin = nb, useWINE=TRUE, OpenBUGS.pgm=OpenBUGS.pgm,
                     WINE=WINE, WINEPATH=WINEPATH,
                     debug = TRUE) 

print(out_siteyr_r, 3)