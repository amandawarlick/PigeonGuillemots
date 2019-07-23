library(dplyr)
library(rvest)
library(readr)
library(lubridate)
library(zoo)
library(NOAATides)
library(ggplot2)

setwd('~/Documents/SAFS/PigeonGuillemots')

#install_github("warlicks/NOAATides")

# NOAA Tides and temps API; https://tidesandcurrents.noaa.gov/tide_predictions.html

# Function to grab from API; https://tidesandcurrents.noaa.gov/api/
# query_tides_data <- function(station_id,
#                              start_date,
#                              end_date,
#                              data_product,
#                              units = 'english',
#                              time_zone = 'lst_ldt',
#                              datum = NULL,
#                              interval = NULL,
#                              verbose = FALSE) {
#   base_url <- "https://tidesandcurrents.noaa.gov/api/datagetter"
#   
#   ## Create a list of all params
#   query_params <- list(station = station_id,
#                        begin_date = start_date,
#                        end_date = end_date,
#                        product = data_product,
#                        datum = datum,
#                        units = units,
#                        time_zone = time_zone,
#                        interval = interval,
#                        format = 'json')
#   # Set up the full url
#   query_url <- httr::modify_url(base_url, query = query_params)
#   if (verbose) {
#     print(query_url)
#   }
#   api_call <- httr::GET(query_url)
#   parsed <- httr::content(api_call, as = 'text')
#   df_list <- jsonlite::fromJSON(parsed, simplifyDataFrame = TRUE, flatten = TRUE)
#   df <- df_list[[1]]
# }

#### Temps ######

# id_temp <- c(9444900, 9447130, 9446484)
# name_temp <- c("Port Townsend", "Seattle", 'Tacoma')
# id_temp <- c(9444900, 9447130)
# name_temp <- c("Port Townsend", "Seattle")

#do separately since certain years don't function and then it messes up the loop
#PT
id_temp <- c(9444900)
name_temp <- c("Port Townsend")
years <- c(2008, 2009,
           2010, 
           2011:2018)

all_june <- data.frame()
all_july <- data.frame()
all_aug <- data.frame()
all_sep <- data.frame()

all_years <- data.frame()
for (j in years) {
for (i in id_temp) {
  temps <- query_tides_data(
    station = i,
    start_date = paste(j,'0601 15:00', sep = ''),
    end_date = paste(j, '0630 17:00', sep = ''),
    data_product = 'water_temperature',
    interval = 'h')
  temps$station <- i
  all_june <- rbind(all_june, temps)
  } #june
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j, '0701 15:00', sep = ''),
      end_date = paste(j, '0731 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h') 
    temps$station <- i
    all_july <- rbind(all_july, temps)
  } #july
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j, '0801 15:00', sep = ''),
      end_date = paste(j, '0831 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h')
    temps$station <- i
    all_aug <- rbind(all_aug, temps)
  } #aug
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j, '0901 15:00', sep = ''),
      end_date = paste(j, '0930 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h')
    temps$station <- i
    all_sep <- rbind(all_sep, temps)
  } #sept
  all_years <- rbind(all_june, all_july, all_aug, all_sep) %>% select(-f)
} #years

setwd('~/Documents/SAFS/PigeonGuillemots/WhidbeyData')

sst_pt2010 <- read.csv('CO-OPS_PT2010.csv', header = T, stringsAsFactors = F) %>%
  transform(station = 9444900) %>%
  transform(t = mdy(t))
  
sst_pt <- all_years %>%
  transform(t = as.Date(ymd_hm(t), format = '%Y-%m-%d')) %>%
  transform(v = as.numeric(v)) %>%
  #bind_rows(sst_pt2010) %>%
  rename(date = t) %>%
  transform(month = month(date)) %>%
  transform(day = day(date)) %>%
  transform(year = factor(year(date))) %>%
  transform(yday = yday(date)) %>%
  rename(temp = v) %>%
  filter(temp < 58)

#test <- filter(year == 2010) %>% summarize(cnt = n_distinct(yday))

sst_day_pt <- sst_pt %>%
  group_by(day, station, month, year, yday) %>%
  summarize(temp = mean(temp, na.rm = T))

########seattle
id_temp <- c(9447130)
name_temp <- c("Seattle")
years <- c(2008, 2009, 
           #2010,
           2011:2015, 
           #2016,
           2017:2018)

all_june <- data.frame()
all_july <- data.frame()
all_aug <- data.frame()
all_sep <- data.frame()

all_years_sea <- data.frame()
for (j in years) {
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j,'0601 15:00', sep = ''),
      end_date = paste(j, '0630 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h')
    temps$station <- i
    all_june <- rbind(all_june, temps)
  } #june
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j, '0701 15:00', sep = ''),
      end_date = paste(j, '0731 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h') 
    temps$station <- i
    all_july <- rbind(all_july, temps)
  } #july
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j, '0801 15:00', sep = ''),
      end_date = paste(j, '0831 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h')
    temps$station <- i
    all_aug <- rbind(all_aug, temps)
  } #aug
  for (i in id_temp) {
    temps <- query_tides_data(
      station = i,
      start_date = paste(j, '0901 15:00', sep = ''),
      end_date = paste(j, '0930 17:00', sep = ''),
      data_product = 'water_temperature',
      interval = 'h')
    temps$station <- i
    all_sep <- rbind(all_sep, temps)
  } #sept
  all_years_sea <- rbind(all_june, all_july, all_aug, all_sep) %>% select(-f)
} #years

setwd('~/Documents/SAFS/PigeonGuillemots/PiGuData/EnvData')

sst_sea2010 <- read.csv('CO-OPS_Seattle2010.csv', header = T, stringsAsFactors = F) %>%
  transform(station = 9447130) %>%
  transform(t = mdy(t))
sst_sea2016 <- read.csv('CO-OPS_Seattle2016.csv', header = T, stringsAsFactors = F) %>%
  transform(station = 9447130) %>%
  transform(t = mdy(t))

sst_sea <- all_years_sea %>%
  transform(t = as.Date(ymd_hm(t), format = '%Y-%m-%d')) %>%
  transform(v = as.numeric(v)) %>%
  bind_rows(sst_sea2010) %>%
  bind_rows(sst_sea2016) %>%
  rename(date = t) %>%
  transform(month = month(date)) %>%
  transform(day = day(date)) %>%
  transform(year = factor(year(date))) %>%
  transform(yday = yday(date)) %>%
  rename(temp = v) %>%
  filter(temp < 58)

sst_day_sea <- sst_sea %>%
  group_by(day, station, month, year, yday) %>%
  summarize(temp = mean(temp, na.rm = T))

sst_day <- rbind(sst_day_pt, sst_day_sea)


#### Tides ######

col_id_SS <- data.frame(site = 
                          c("Totten Legacy", "Totten Elizan C", "Totten Elizan A Tower", "Totten Elizan A North",
                              "Youngs Cove A", "Youngs Cove B", "Youngs Cove C", "Youngs Cove D", "Youngs Cove E", "Youngs Cove F", "Flapjack",
                            "Priest Point", "Hearthfire",
                            "Gull Harbor", "Burfoot Park", "Burfoot Extended",
                            "Zangle Cove A", "Big Fish Trap", "Briscoe Pt",
                            "Zittels Marina", "Amsterdam Bay", "Mill Bight A", "Mill Bight B", "Mill Bight Extended", "Walnut Rd A", "Walnut Rd B", 
                              "Walnut Rd C", "Beachcrest", "Butterball Cove South A", "Butterball North", 
                            "Andys Marine Park West", "Andys Marine Park South", "Higgins Cove",
                            "Sandy Pt", "Lyle Point", "Cole Pt",
                            "Ketron Ferry", "Ketron SE", "Ketron SW",
                            "Jarrell Cove State Park", "Harstine Pointe Lagoon", "Compass Rock", "Hope Island NE"), 
                        grp = c(rep(11, 11),
                                12, 12, 13, 13, 13,
                                14, 14, 14, rep(15, 11),
                                rep(16, 3), rep(17, 3),
                                18, 18, 18, rep(19, 4)))
station_id_SS <- data.frame(station = c(9446742, 9446969, 9446807, 9446800, 9446752,
                                        9446671, 9446804, 9446714, 9446489), 
                            grp = 11:19) %>%
  merge(col_id_SS, by = 'grp') %>%
  transform(region = 'SS')

id_SS <- c(9446742, 9446969, 9446807, 9446800, 9446752,
           9446671, 9446804, 9446714, 9446489)

## Whidbey and South Sound station numbers and corresponding colony names

col_id_whid <- data.frame(site = 
                       c("Crescent Harbor", "Forbes Point", "Maylor Point",
                         "Coupeville Wharf", "Monroe Landing", "Rolling Hills",
                         "Harrington North", "Harrington South", "Pratts Bluff",
                         "Langley Marina",
                         "Possession Point",
                         "Bush Point Dock", "Malmo Bluff", "Shore Meadows", "Mutiny Sands", "Limpet Lane", 
                            "Double Bluff North", "Double Bluff South", 
                         "Lagoon North #0", "Lagoon North #1", "Lagoon North #2", "Lagoon North #3", "Lagoon South", 
                            "Keystone", "Ledgewood", "Hancock Lake", "Fort Casey", 
                         "Libby North",
                         "Swantown", "Hastie Lake", "Hastie Lake North",
                         "Cliffside"), 
                     grp = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, rep(6, 7),
                             rep(7, 9),
                             8, 9, 9, 9, 10))
station_id_whid <- data.frame(station = c(9447952, 9447929, 9447883, 9447856, 9447814, 9447854, 9447905, 
                                     9447934, 9447951, 9447973), grp = 1:10) %>%
  merge(col_id_whid, by = 'grp') %>%
  transform(region = 'Whidbey') 

id_whid <- c(9447952, 9447929, 9447883, 9447856, 9447814, 9447854, 9447905, 
             9447934, 9447951, 9447973)

all_whidbey <- data.frame()

for (i in id_whid) {
  tides <- query_tides_data(
    station = i,
    start_date = '20080601',
    end_date = '20150930',
    data_product = 'predictions',
    datum = 'MLLW',
    interval = 'hilo')
  tides$station <- i
  all_whidbey <- rbind(all_whidbey, tides)
}

all_whidbey2 <- data.frame() #break into 2 so don't exceed query quota

for (i in id_whid) {
  tides <- query_tides_data(
    station = i,
    start_date = '20160601',
    end_date = '20180930',
    data_product = 'predictions',
    datum = 'MLLW',
    interval = 'hilo')
  tides$station <- i
  all_whidbey2 <- rbind(all_whidbey2, tides)
}

all_whidbey <- all_whidbey %>%
  bind_rows(all_whidbey2) %>%
  merge(station_id_whid %>% select(-grp), by = 'station')

all_SS <- data.frame()

for (i in id_SS) {
  tides <- query_tides_data(
    station = i,
    start_date = '20150601',
    end_date = '20180930',
    data_product = 'predictions',
    datum = 'MLLW',
    interval = 'hilo')
  tides$station <- i
  all_SS <- rbind(all_SS, tides)
}

all_SS <- all_SS %>%
  merge(station_id_SS %>% select(-grp), by = 'station')

# tides_all <- all_SS %>%
#   bind_rows(all_whidbey) %>%
#   transform(month = month(t), date = as.Date(t)) %>%
#   filter(month >=6 & month <= 9) %>%
#   transform(time = as.POSIXlt(t, "%Y-%m-%d %H:%M:%S"))

tides_all <- all_whidbey %>%
  transform(month = month(t), date = as.Date(t)) %>%
  filter(month >=6 & month <= 9) %>%
  transform(time = as.POSIXlt(t, "%Y-%m-%d %H:%M:%S"))

setwd('~/Documents/SAFS/PigeonGuillemots')

write.csv(tides_all, "tides_all.csv", row.names = F)

write.csv(sst_day, 'sst_day.csv', row.names = F)


