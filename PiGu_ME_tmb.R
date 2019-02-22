
library(TMB)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("~/Documents/SAFS/PigeonGuillemots")

set.seed(1)

#####################################

# read in and format data
data_MEJS <- read.csv("burrow_CH.csv", header = T, stringsAsFactors = F)
#data_MEJS <-subset(data_MEJS, site == "Double Bluff North" &  year == 2016)
#data_MEJS <-subset(data_MEJS, site == "Double Bluff North")

data_MEJS_wide <- data_MEJS %>% 
  distinct(region, year, site, burrow_name, study_day, capt_hist) %>% 
  dcast(region + year + site + burrow_name ~ study_day, value.var = 'capt_hist', fun.aggregate = mean)

CH_MEJS <- as.matrix(data_MEJS_wide %>% select(-c(region, year, site, burrow_name)), rownames.force = F)
colnames(CH_MEJS) <- NULL
CH_MEJS <- round(CH_MEJS, 0)

# Augment data
# nz <- 500
# CH.aug <- rbind(CH_MEJS, matrix(3, nrow = nz, ncol = dim(CH_MEJS)[2])) #3 for not seen
# 
# #augmented effort - not finished
# eff.aug <- rbind(EFFORT_MEJS, matrix(EFFORT_MEJS[1,], nrow = nz, ncol = dim(EFFORT_MEJS)[2], byrow=T))
# eff.aug <- eff.aug[,-1] #remove first column (dummy column), because detection goes off of t-1?

## Create binary effort matrix with 1 if site surveyed on day t, 0 otherwise.
eff <- matrix(0, nr = dim(CH_MEJS)[1], nc = dim(CH_MEJS)[2])
eff[!is.na(CH_MEJS)] <- 1

# define various quantities
Nind <- dim(CH_MEJS)[1]
N_occ <- dim(CH_MEJS)[2]

######################################

get.first <- function(x) min(which(x < 3))
f <- apply(CH_MEJS, 1, get.first)
#a <- match(Inf, f)

get.last <- function(x) max(which(x < 3))
l <- apply(CH_MEJS, 1, get.last)


# transpose data
#data <- t(data)

###############################################

compile("PiGu_ME_tmb.cpp")
dyn.load(dynlib("PiGu_ME_tmb"))

model <- MakeADFun(
  data = list(ch = data, fc = fc, fs = init.state), 
  parameters = list(),
  DLL = "PiGu_ME_tmb")
opt <- do.call("optim", model) # optimisation
model$fn(binit) # evaluate likelihood at the inits
model$report()$B # display B
model$report()$BE # display BE
model$report()$A # display A
model$report()$PROP # display PROP
rep <- sdreport(model)
rep # get SEs


####################################################

# #same model in R code
# devMULTIEVENT <- function(b,data,eff,e,garb,nh,km1){
# 
# # data encounter histories, eff counts
# # e vector of dates of first captures
# # garb vector of initial states 
# # km1 nb of recapture occasions (nb of capture occ - 1)
# # nh nb ind
# 
# # OBSERVATIONS (+1)
# # 0 = non-detected
# # 1 = seen and ascertained as non-breeder
# # 2 = seen and ascertained as breeder
# # 3 = not ascertained
# 
# # STATES
# # 1 = alive non-breeder
# # 2 = alive breeder
# # 3 = dead
# 
# # PARAMETERS
# # phiNB  survival prob. of non-breeders
# # phiB  survival prob. of breeders
# # pNB  detection prob. of non-breeders
# # pB  detection prob. of breeders
# # psiNBB transition prob. from non-breeder to breeder
# # psiBNB transition prob. from breeder to non-breeder
# # piNB prob. of being in initial state non-breeder
# # deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
# # deltaB prob to ascertain the breeding status of an individual encountered as breeder
# 
# # logit link for all parameters
# # note: below, we decompose the state and obs process in two steps composed of binomial events, 
# # which makes the use of the logit link appealing; 
# # if not, a multinomial (aka generalised) logit link should be used
# par = plogis(b)
# piNB <- par[1]
# phiNB <- par[2]
# phiB <- par[3]
# psiNBB <- par[4]
# psiBNB <- par[5]
# pNB <- par[6]
# pB <- par[7]
# deltaNB <- par[8]
# deltaB <- par[9]
# 
# # prob of obs (rows) cond on states (col)
# B1 = matrix(c(1-pNB,pNB,0,1-pB,0,pB,1,0,0),nrow=3,ncol=3,byrow=T)
# B2 = matrix(c(1,0,0,0,0,deltaNB,0,1-deltaNB,0,0,deltaB,1-deltaB),nrow=3,ncol=4,byrow=T)
# B = t(B1 %*% B2)
# 
# # first encounter
# BE1 = matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=T)
# BE2 = matrix(c(1,0,0,0,0,deltaNB,0,1-deltaNB,0,0,deltaB,1-deltaB),nrow=3,ncol=4,byrow=T)
# BE = t(BE1 %*% BE2) 
# 
# # prob of states at t+1 given states at t
# A1 <- matrix(c(phiNB,0,1-phiNB,0,phiB,1-phiB,0,0,1),nrow=3,ncol=3,byrow=T)
# A2 <- matrix(c(1-psiNBB,psiNBB,0,psiBNB,1-psiBNB,0,0,0,1),nrow=3,ncol=3,byrow=T)
# A <- A1 %*% A2
# 
# # init states
# PI <- c(piNB,1-piNB,0)
# 
# # likelihood
# l <- 0
# for (i in 1:nh) # loop on ind
# {
#   ei <- e[i] # date of first det
#   oe <- garb[i] + 1 # init obs
#   evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
#   ALPHA <- PI*BE[oe,]
#   for (j in (ei+1):(km1+1)) # cond on first capture
#   {
#   if ((ei+1)>(km1+1)) {break}
#   ALPHA <- (ALPHA %*% A)*B[evennt[j],]
#   }
#   l <- l + log(sum(ALPHA))#*eff[i]
# }
# l <- -l
# l
# }
# 
# #benchmarking
# library(microbenchmark)
# res = microbenchmark(
#   optim(binit,devMULTIEVENT,NULL,hessian=F,data,eff,fc,init.state,nh,km1,method="BFGS"),
#   do.call("optim", f),
#   times=5
# )
# res2 = summary(res)

  
  