
library(tidyr)
#library(ggfortify)
library(broom) #tidy() augment()
#library(lme4)
#library(lmerTest)
library(scales)
library(reshape2)
library(stats) 
library(sciplot)
library(jagsUI)
library(chron)
library(rstudioapi)
library(rjags)
library(jagsUI)
library(R2OpenBUGS)
library(dplyr)
library(data.table)
#library(ggplot2)
library(loo)
library(knitr)
library(stringr)
library(lubridate)
#library(IDPmisc)
library(zoo)
library(magrittr)


setwd("~/Documents/SAFS/PigeonGuillemots")
path <- getActiveDocumentContext()$path 
setwd(dirname(path))

#contains single row where burrow_name == NA when no activity was detected, but survey started
data_burrow_all <- read.csv("data_burrow.csv", header = T, stringsAsFactors = F) %>%
  filter(region == 'Whidbey') %>%
  dplyr::select(-c(Gunnel, Sculpin, Other)) %>%
  filter(site %in% c('Mutiny Sands', 'Double Bluff North')) %>%
  filter(year > 2016) %>%
  distinct()  #2 duplicates, all lagoon south 2017

#study start day per year, island-wide
start_end_all <- data_burrow_all %>%
  group_by(year, region) %>%
  summarize(start_day = min(yday, na.rm = T), end_day = max(yday, na.rm = T))

no_vis <- data_burrow_all %>%
  filter(is.na(burrow_name)) %>%
  dplyr::select(region, year, site, week, yday)

all_days <- data_burrow_all %>%
  #filter(is.na(burrow_name)) %>%
  group_by(region, year, site) %>%
  distinct(yday)
fill <- data_burrow_all %>%  #all day-burrow combinations
  filter(!is.na(burrow_name)) %>%
  group_by(region, year, site, burrow_name) %>%
  distinct(burrow_name) %>%
  merge(all_days, all = T) 

burrow <- data_burrow_all %>%
  filter(!is.na(burrow_name)) %>% #remove all dummy rows where no burrows were seen
  merge(fill, by = c('region', 'site', 'year', 'burrow_name', 'yday'), all = T) %>%
  merge(start_end_all, by = c('region', 'year'), all = T) %>%
  transform(study_day = yday - start_day + 1) #%>% arrange(burrow_name)

#burrow visit and prey visit start stops
prey_start_end <- burrow %>%
  filter(tot_prey > 0) %>%
  group_by(region, year, site, burrow_name) %>%
  summarize(prey_start = min(study_day), prey_end = max(study_day))  

bv_start_end <- burrow %>%
  filter(burrow_visit > 0) %>%
  group_by(region, year, site, burrow_name) %>%
  summarize(bv_start = min(study_day), bv_end = max(study_day))  

start_end_visits <- prey_start_end %>%
  merge(bv_start_end, by = c('region', 'year', 'site', 'burrow_name'), all = T)

day_range <- burrow %>%
  group_by(region, year, site) %>%
  summarize(min_day = min(study_day, na.rm = T), max_day = max(study_day, na.rm = T))

week_range <- burrow %>%
  group_by(region, year, site) %>%
  summarize(min_week = min(week), max_week = max(week))

n_visits <- burrow %>%
  group_by(region, year, site) %>% #don't include burrow here, since hasn't been expanded yet
  summarize(n_visits = n_distinct(study_day)) #use study_day instead of date - NAs in date

#combining all data with dummy framework of all days, create capture history
burrow_CH <- burrow %>%
  merge(start_end_visits, by = c('region', 'year', 'site', 'burrow_name'), all.x = T) %>%
  transform(prey_days = prey_end + 1 - prey_start) %>%
  transform(bv_days = bv_end + 1 - bv_start) %>%
  group_by(region, year, site) %>%
  merge(day_range, by = c('region', 'year', 'site')) %>%
  merge(n_visits, by = c('region', 'year', 'site')) %>% #arrange(burrow_name)
  transform(capt_hist = ifelse(is.na(burrow_visit) & is.na(tot_prey), 3, #observed but not detected
                               ifelse(burrow_visit == 0 & tot_prey == 0, 3, #observed but not detected
                                      ifelse(tot_prey > 0, 2, #prey visit
                                             ifelse(burrow_visit > 0, 1, #burrow visit
                                                    100))))) %>% 
  select(region, year, site, week, yday, start_day, study_day, n_visits, 
         min_day, max_day, burrow_name, capt_hist) %>% distinct() %>%
  dcast(region + year + site + burrow_name + min_day + max_day ~ study_day, value.var = 'capt_hist', 
        fun.aggregate = mean) %>%
  filter(!is.na(burrow_name)) %>%
  melt(id.vars = c('region', 'year', 'site', 'burrow_name', 'min_day', 'max_day')) %>%
  transform(study_day = as.numeric(as.character(variable))) %>% 
  transform(capt_hist = ifelse(study_day < min_day | study_day > max_day, NA, value)) %>%
  dplyr::select(region, year, site, burrow_name, study_day, capt_hist) %>% #NA outside start/stop days
  dcast(region + year + site + burrow_name ~ study_day, value.var = 'capt_hist') 

#df of survey days at burrow-week level
survey <- burrow %>%
  merge(n_visits, by = c('region', 'year', 'site')) %>%
  merge(week_range, by = c('region', 'year', 'site')) %>%
  #transform(occ = week-min_week) %>%
  dplyr::select(region, year, site, yday, burrow_name, study_day, n_visits, week, max_week) %>%
  #dcast(region + year + site + burrow_name + n_visits ~ week, value.var = 'study_day', fun.aggregate = mean)
  dcast(region + year + site + n_visits + burrow_name ~ study_day, value.var = 'study_day', fun.aggregate = mean)

max_recap <- max(survey$n_visits) #max number of resights for ncol, NOT distinct study days

####augment

#number of nests per site per year - for proportion
num_nests <- burrow %>%
  group_by(region, year, site) %>%
  summarize(nest_cnt = n_distinct(burrow_name)) %>%
  transform(burrow_name = 'XXX') %>%
  merge(n_visits, by = c('region', 'year', 'site')) %>%
  transform(aug = round(0.2*nest_cnt)) %>%
  transform(aug = ifelse(aug < 1, 1, aug)) %>%
  dplyr::select(-nest_cnt)

aug_fun <- function(df, id_cols, visit_col){
  # Keep only the needed columns. 
  cols <- append(id_cols, visit_col)
  df_down_select <- df %>% select(cols)
  # Create a storage data frame
  final_df <- data.frame(stringsAsFactors = FALSE)
  # Loop through the data frame
  for (i in 1:nrow(df_down_select)){
    # Replicate the id data based on the number of visits
    list_rep <- rep(df_down_select[i, id_cols], df_down_select[i, visit_col])
    # Convert the resulting list to a data frame
    current_df <- as.data.frame(matrix(unlist(list_rep),
                                       nrow = df_down_select[i, visit_col],
                                       byrow = TRUE),
                                stringsAsFactors = FALSE)
    final_df <- bind_rows(final_df, current_df)
  }
  return(final_df)
}

aug_id <- aug_fun(num_nests, c('region', 'year', 'site', 'burrow_name', 'n_visits'), 'aug') %>%
  transform(V4 = 'XXX')
colnames(aug_id) <- c('region', 'year', 'site', 'burrow_name', 'n_visits') 

#take just id cols and days
days_obs <- survey %>%
  dplyr::select(-c(burrow_name)) %>%
  distinct()

#merge id columns and day observation columns
aug_inds <- aug_id %>%
  merge(days_obs, by = c('region', 'year', 'site', 'n_visits'), all = T) %>%
  transform(year = as.numeric(year))
colnames(aug_inds) <- colnames(survey)

#take augmented survey day data and turn to '3' obs
ch_obs <- days_obs
for (i in 1:nrow(ch_obs)) {
  for (j in 5:ncol(ch_obs)) {
    if (!is.nan(ch_obs[i,j])) {
      ch_obs[i,j] <- 3
    }
  }
}

aug_inds_CH <- aug_id %>%
  merge(ch_obs, by = c('region', 'year', 'site', 'n_visits'), all = T) %>%
  transform(year = as.integer(year)) %>%
  dplyr::select(-n_visits)
colnames(aug_inds_CH) <- colnames(burrow_CH)

burrow_CH_aug <- burrow_CH %>%
  bind_rows(aug_inds_CH) 

##collapse so that all data are at beginning, trailing NAs for the rest
ch_dat <- matrix(NA, nrow = dim(burrow_CH_aug)[1], ncol = max_recap)

for (i in 1:dim(burrow_CH_aug)[1]) {
  temp <- as.numeric(burrow_CH_aug[i, 4 + c(which(!is.na(burrow_CH_aug[i,5:dim(burrow_CH_aug)[2]])))])
  num <- length(temp)
  ch_dat[i,] <- c(temp, rep(NA, max_recap-num))
}

#for unique year-site combination, need vector of study_day that is length of n_visits

max_visits <- max(n_visits$n_visits)

surv_aug <- survey %>%
  bind_rows(aug_inds %>% transform(n_visits = as.numeric(n_visits)))

survey <- surv_aug

survey_days <- matrix(NA, nrow = dim(survey)[1], ncol = max_visits)

for (i in 1:dim(survey)[1]) {
  temp <- as.numeric(survey[i, 5 + c(which(!is.na(survey[i,6:dim(survey)[2]])))])
  #num <- survey[i,5]
  num <- length(temp)
  survey_days[i,] <- c(temp, rep(NA, max_visits-num))
}

year_i <- as.numeric(factor(burrow_CH_aug$year))
site_i <- survey$site
visits_i <- survey$n_visits

#awkward renaming of things from different scripts
effort <- survey_days
#ch <- read.csv('ch_dat.csv', header = T)
ch <- ch_dat
nyear <- length(unique(year_i))
y <- ch

#model
# -------------------------------------------------
# Parameters:
# phiA: survival probability from egg to chick
# phiB: survival probability from chick to fledge
# psiAB: probability of transitioning from egg to chick
# pA: detection probability of egg burrow
# pB: detection probability of chick burrow
# b: conditional on there being a chick, probability of seeing just a burrow visit
# gamma: entry probability
# alpha: conditional on entry at occasion 1, probability burrow had a chick. alpha set to 0 for t>1. So if a burrow is initiated after day one, it must start as an egg burrow.  

# -------------------------------------------------
# States:
# 1 not entered
# 2 alive as egg
# 3 alive as chick
# 4 terminated; dead or fledged

# Observations:  
# 1 Burrow visit
# 2 Prey visit
# 3 not seen
# NA not observed - no survey that day
# -------------------------------------------------


model_MEJS <- function () {
  
  # Priors and constraints
  for (i in 1:M) {
    for (t in 1:(n.occasions-1)){
      logit(phiA[i,t]) <- mu.phi[1] + b.egg.y[year_i[i]] 
      logit(phiB[i,t]) <- mu.phi[2] + b.chick.y[year_i[i]] 
      gamma[i,t] ~ dunif(0, 1)    # Prior for entry probabilities at occasion t
      pA[i,t] <- mean.p[1]     # egg burrow detection
      pB[i,t] <- mean.p[2]     # chick burrow detection
      psiAB[i,t] <- mean.psiAB
    } #t
  } #i
  
  b ~ dunif(0,1)            # prior for assignment probability
  mean.psiAB ~ dunif(0,1)   # transition from egg to chick
  pi ~ ddirch(alpha[1:3])   #dirichlet prior for multinomial
  alpha[1] <- 1/3
  alpha[2] <- 1/3
  alpha[3] <- 1/3
  
  for (u in 1:2){
    mu.phi[u] ~ dnorm(0, 0.001) 
    #mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
    mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
  }
  
  for (y in 2:nyear) {
    b.chick.y[y] ~ dnorm(0, 0.001)
    b.egg.y[y] ~ dnorm(0, 0.001)
  }
  b.chick.y[1] <- 0
  b.egg.y[1] <- 0
  
  #inside year loop
  phiA.int <- 1 / (1+exp(-mu.phi[1]+beta.year[y]))
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
      ps[1,i,t,1] <- 1-gamma[i,t]                    #probability of not entering
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
  for (t in 2:n.occasions) {
    chick[t] <- sum(z[1:M,t] == 3)
    egg[t] <- sum(z[1:M,t] == 2)
  } #t
  
##annual level, help
  for (y in year_i) {
    for (i in 1:M) {
      days.chick.y[i,y] <- sum(z[i,] == 3)
    fledged.low.y[i,y] <- step(days.chick.y[i]-45) 
    } #i
    #   everActive.y[y] <- max(z[1:M,]>1) #need the number of active ever
    n.fledged.y[y] <- sum(fledged.low.y[1:M,y])
} #y

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

# Bundle data
jags.data <- list(y = y, n.occasions = max(effort, na.rm = T), nyear = length(unique(year_i)),
                  #z = z.st, 
                  effort = effort, #from setup, survey_days
                  visits_i = visits_i, year_i = year_i,
                  M = dim(ch)[1])

inits <- function(){list(#mean.phi = runif(2, 0, 1), 
  #z = z.init,
  z = matrix(3, nrow = dim(ch)[1], ncol = max(effort, na.rm = T)),
  mean.p = runif(2, 0, 1))}  

# Parameters monitored
parameters <- c('phiA.int', 'phiB.int', "mean.p", "b", 'fledged.low.y',
                'n.active.burrow', 'n.fledged.low', 'n.fledged.high', 'mean.days.chick', 'mean.days.egg',
                #"gamma", 'z', 'fledged', 
                'b.egg.y', 'b.chick.y',
                "mean.psiAB", "N.egg", "N.chick", 'nest.succ.low', 'nest.succ.high', 
                "N.active", "Nstar", "psi", 'chick', 'egg', 'gamma')

# MCMC settings
ni <- 8; nt <- 1; nb <- 5; nc <- 2

out_pigu <- jags(jags.data, inits, parameters, model.file = model.file,
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)



