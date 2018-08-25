library(dplyr)
library(rvest)
library(readr)
library(lubridate)
library(zoo)

# station IDs


#############################
# NOAA Tides; https://tidesandcurrents.noaa.gov/tide_predictions.html

# Function to grab from API; https://tidesandcurrents.noaa.gov/api/
query_tides_data <- function(station_id,
                             start_date,
                             end_date,
                             data_product,
                             units = 'english',
                             time_zone = 'lst_ldt',
                             datum = NULL,
                             interval = NULL,
                             verbose = FALSE){
  base_url <- "https://tidesandcurrents.noaa.gov/api/datagetter"
  
  ## Create a list of all params.
  query_params <- list(station = station_id,
                       begin_date = start_date,
                       end_date = end_date,
                       product = data_product,
                       datum = datum,
                       units = units,
                       time_zone = time_zone,
                       interval = interval,
                       format = 'json')
  # Set up the full url
  query_url <- httr::modify_url(base_url, query = query_params)
  if (verbose) {
    print(query_url)
  }
  api_call <- httr::GET(query_url)
  parsed <- httr::content(api_call, as = 'text')
  df_list <- jsonlite::fromJSON(parsed, simplifyDataFrame = TRUE, flatten = TRUE)
  df <- df_list[[1]]
}


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
    end_date = '20170930',
    data_product = 'predictions',
    datum = 'MLLW',
    interval = 'hilo')
  tides$station <- i
  all_whidbey <- rbind(all_whidbey, tides)
}

all_whidbey <- all_whidbey %>%
  merge(station_id_whid %>% select(-grp), by = 'station')

all_SS <- data.frame()

for (i in id_SS) {
  tides <- query_tides_data(
    station = i,
    start_date = '20150601',
    end_date = '20170930',
    data_product = 'predictions',
    datum = 'MLLW',
    interval = 'hilo')
  tides$station <- i
  all_SS <- rbind(all_SS, tides)
}

all_SS <- all_SS %>%
  merge(station_id_SS %>% select(-grp), by = 'station')

tides_all <- all_SS %>%
  bind_rows(all_whidbey) %>%
  transform(month = month(t), date = as.Date(t)) %>%
  filter(month >=6 & month <= 9) %>%
  transform(time = as.POSIXlt(t, "%Y-%m-%d %H:%M:%S"))

write.csv(tides_all, "tides_all.csv", row.names = F)


######## Archived Stuff

####### #Tides.net for 2014-2018 fully functional, keep to compare

# ## Whidbey
# # 644 <- CrescentHarbor (MarinersCove, CrescentHarbor, ForbesPoint, MaylorPoint)
# # 620 <- Coupeville (MonroeLanding, RollingHills, Coupeville)
# # 1078 <- Greenbank (HarringtonN, HarringtonS, PrattsBluff)
# # 2395 <- SandyPoint (Langley)
# # 1026 <- Glendale (PossessionPoint)
# # 361 <- BushPoint (Bush Point Dock, MalmoBluff, ShoreMeadows, MutinySands, LimpetLane, DoubleBluff
# #                   LagoonN, LagoonS, Keystone, Ledgewood, Hancock, FortCasey)
# # 2686 <- SunsetBeach (Cliffside, Swantown, HastieLake)
# 
# col_id <- data.frame(site = 
#         c("Crescent Harbor", "Forbes Point", "Maylor Point",
#           "Coupeville Wharf", "Monroe Landing", "Rolling Hills",
#           "Harrington North", "Harrington South", "Pratts Bluff",
#           "Langley Marina",
#           "Possession Point",
#           "Bush Point Dock", "Malmo Bluff", "Shore Meadows", "Mutiny Sands", "Limpet Lane", 
#                 "Double Bluff North", "Double Bluff South", "Lagoon North #0", "Lagoon North #1", 
#                 "Lagoon North #2", "Lagoon North #3", "Lagoon South", 
#                 "Keystone", "Ledgewood", "Hancock Lake", "Fort Casey", "Libby North",
#           "Cliffside", "Swantown", "Hastie Lake", "Hastie Lake North"), 
#                              grp = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, rep(6, 17),
#                                      7, 7, 7, 7))
# station_id <- data.frame(station = c(644, 620, 1078, 2395, 1026, 361, 2686), grp = 1:7) %>%
#   merge(col_id, by = 'grp')

# id <- c(644, 620, 1078, 2395, 1026, 361, 2686)
# month <- formatC(5:9, width = 2, flag = 0)
# year <- 2014:2018
# 
# url_all = data.frame()
# 
# for (k in id) {
#   for (i in year) {
#     url <- paste0("https://www.tides.net/washington/", k, "/?year=",
#                   i, "&month=", month)
#     #print(url)
#     df <- data.frame(url, stringsAsFactors = F)
#     url_all <- rbind(url_all, df)
#     }
# }
# 
# urls <- url_all[,1]

#read html using rvest()
#tides_all <- data.frame()

# for (i in urls) {
# tides <- read_html(i) %>%
#   html_nodes("table") %>%
#   html_table(fill = T)
#   if (nrow(tides[[35]]) > 5) {
#       tides <- tides[[35]]
#         } else { 
#       tides <- tides[[34]]
#     }
# tides$url <- i
# colnames(tides) <- c('date', rep(c('dif', 'type', 'tide', 'time'), 4), 'dif', 'sun', 'moon', 'url')
# #print(tides)
#  df <- data.frame(tides)
#  tides_all <- rbind(tides_all, df)
# }
# 
# tides_all_clean <- tides_all %>%
#   transform(station = sub(".*n/", "", url)) %>%
#   transform(station = sub("/.*", "", station)) %>%
#   transform(year = sub(".*r=", "", url)) %>%
#   transform(year = sub("&m.*", "", year)) %>%
#   transform(date = as.Date(paste(date, year, sep = "/"), format = "%m/%d/%Y")) %>%
#   transform(time = sub(".* ", "" ,strptime(time, format = "%I:%M%p")),
#             time.1 = sub(".* ", "" ,strptime(time.1, format = "%I:%M%p"))) %>%
#   transform(time = as.POSIXlt(paste(date, time, sep = " ")),
#             time.1 = as.POSIXlt(paste(date, time.1, sep = " "))) %>%
#   transform(tide_ft = as.numeric(sub("'.*", "", tide)),
#             tide_in = as.numeric(parse_number(sub(".*'", "", tide)))/12,
#             tide_ft.1 = as.numeric(sub("'.*", "", tide.1)),
#             tide_in.1 = as.numeric(parse_number(sub(".*'", "", tide.1)))/12) %>%
#   transform(tide = tide_ft + tide_in,
#             tide.1 = tide_ft.1 + tide_in.1) %>%
#   select(-matches("_"))
# 
# tides_tidy <- tides_all_clean %>%
#   select(year, station, date, type, tide, time, type.1, tide.1, time.1) %>%
#   transform(sample = as.POSIXlt(paste(date, 
#               format(strptime("8:00", format = "%H:%M"), "%H:%M:%S"), sep = " "))) %>%
  # transform(from_low =
  #             ifelse(type == "L", as.numeric(time-sample, units = "mins"),
  #                    as.numeric(time.1-sample, units = "mins")),
  #           from_high =
  #             ifelse(type == "H", as.numeric(time-sample, units = "mins"),
  #                    as.numeric(time.1-sample, units = "mins"))) %>%
  # merge(station_id %>% select(-grp), by = 'station')
#write.csv(tides_tidy, "tides_tidy.csv", row.names = F)




###### Messy other things I tried for NOAA tides and threw out

#binding test
# for (i in id) {
#   all_stations <- data.frame()
#   all_one <- data.frame()
#   for (k in b_dates) for (j in e_dates) {
#     all_oneID <- paste0(k, "to", j)
#     print(all_oneID) }
#     all_one <- rbind(all_one, all_oneID)
#   }
#   all_stations <- rbind(all_stations, all)
# }

# bmo <- c("0601", "0701", "0801", "0901")
# emo <- c("0630", "0731", "0831", "0930")
# bmo <- c("0601")
# emo <- c("0930")
# year <- 2008:2017

# url_old = data.frame()
# 
# for (k in id) {
#     url <- paste0("https://tidesandcurrents.noaa.gov/noaatidepredictions.html?id=", k, 
#                   "&units=standard&bdate=", b_dates, 
#     "&edate=", e_dates, "&timezone=LST/LDT&clock=12hour&datum=MLLW&interval=hilo&action=data")
#     df <- data.frame(url, stringsAsFactors = F)
#     url_old <- rbind(url_old, df)
# }

# b_dates = data.frame()
# for (i in year) {
#   bdate <- paste0(i, bmo)
#   df <- data.frame(bdate, stringsAsFactors = F)
#   b_dates <- rbind(b_dates, df)
# }
# b_dates <- as.character(b_dates[,])
# e_dates = data.frame()
# for (i in year) {
#   edate <- paste0(i, emo)
#   df <- data.frame(edate, stringsAsFactors = F)
#   e_dates <- rbind(e_dates, df)
# }
# e_dates <- as.character(e_dates[,])
