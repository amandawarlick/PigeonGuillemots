library(dplyr)
library(rvest)

#station IDs
# 604 <- DeceptionPass (DeceptionPass)
#* 644 <- CrescentHarbor (MarinersCove, CrescentHarbor, ForbesPoint, MaylorPoint)
#* 620 <- Coupeville (MonroeLanding, RollingHills, Coupeville)
#* 1078 <- Greenbank (HarringtonN, HarringtonS, PrattsBluff)
#* 2395 <- SandyPoint (Langley)
#* 1026 <- Glendale (PossessionPoint)
#* 361 <- BushPoint (MalmoBluff, ShoreMeadows, MutinySands, LimpetLane, DoubleBluff
#                   LagoonN, LagoonS, Keystone, Ledgewood, Hancock, FortCasey)
# 16 <- AdmiraltyHead (Keystone, Ledgewood, Hancock, FortCasey)
#* 2686 <- SunsetBeach (Cliffside, Swantown, HastieLake)

col_id <- data.frame(col = c("Deception Pass",
                                     "Mariners Cove", "Crescent Harbor", "Forbes Point", "Maylor Point",
                                     "Coupeville", "Monroe Landing", "Rolling Hills #1", "Rolling Hills #2", "Rolling Hills #3", "Rolling Hills",
                                     "Harrington North", "Harrington South", "Pratts Bluff",
                                     "Langley",
                                     "Possession Point",
                                     "Malmo Bluff", "Shore Meadows", "Mutiny Sands", "Limpet Lane", "Double Bluff", "Lagoon North", "Lagoon South",
                                     "Keystone", "Ledgewood", "Hancock Lake", "Fort Casey",
                                     "Cliffside", "Swantown", "Hastie Lake"), 
                             grp = c(2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 
                                     5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8))

station_id <- data.frame(id = c(604, 644, 620, 1078, 2395, 1026, 361, 16, 2686), grp = 1:9) %>%
  merge(col_id, by = 'grp')

id <- c(644, 620, 1078, 2395, 1026, 361, 2686)
#id <- c(644, 620)
month <- formatC(5:9, width = 2, flag = 0)
year <- 2014:2018

url_all = data.frame()

for (k in id) {
  for (i in year) {
    url <- paste0("https://www.tides.net/washington/", k, "/?year=",
                  i, "&month=", month)
    #print(url)
    df <- data.frame(url, stringsAsFactors = F)
    url_all <- rbind(url_all, df)
    }
}

urls <- url_all[,1]

#read html using rvest()
tides_all <- data.frame()
#url_test <- urls[1:3]

for (i in urls) {
tides <- read_html(i) %>%
  html_nodes("table") %>%
  html_table(fill = T)
  if (nrow(tides[[35]]) > 5) {
      tides <- tides[[35]]
        } else { 
      tides <- tides[[34]]
    }
tides$url <- i
colnames(tides) <- c('date', rep(c('dif', 'type', 'tide', 'time'), 4), 'dif', 'sun', 'moon', 'url')
#print(tides)
 df <- data.frame(tides)
 tides_all <- rbind(tides_all, df)
}

tides_all_clean <- tides_all %>%
  transform(station = sub(".*n/", "", url)) %>%
  transform(station = sub("/.*", "", station)) %>%
  transform(year = sub(".*r=", "", url)) %>%
  transform(year = sub("&m.*", "", year)) %>%
  transform(date = as.Date(paste(date, year, sep = "/"), format = "%m/%d/%Y")) %>%
  transform(time = sub(".* ", "" ,strptime(time, format = "%I:%M%p")),
            time.1 = sub(".* ", "" ,strptime(time.1, format = "%I:%M%p"))) %>% 
  transform(time = as.POSIXlt(paste(date, time, sep = " ")),
            time.1 = as.POSIXlt(paste(date, time.1, sep = " "))) %>%
  transform(tide_ft = as.numeric(sub("'.*", "", tide)),
            tide_in = as.numeric(parse_number(sub(".*'", "", tide)))/12,
            tide_ft.1 = as.numeric(sub("'.*", "", tide.1)),
            tide_in.1 = as.numeric(parse_number(sub(".*'", "", tide.1)))/12) %>%
  transform(tide = tide_ft + tide_in,
            tide.1 = tide_ft.1 + tide_in.1) %>%
  select(-matches("_"))

tides_tidy <- tides_all_clean %>%
  select(year, station, date, type, tide, time, type.1, tide.1, time.1) %>%
  transform(sample = as.POSIXlt(paste(date, format(strptime("8:00", format = "%H:%M"), "%H:%M:%S"), sep = " "))) %>%
  transform(from_low = 
              ifelse(type == "L", as.numeric(time-sample, units = "mins"), as.numeric(time.1-sample, units = "mins")),
            from_high = 
              ifelse(type == "H", as.numeric(time-sample, units = "mins"), as.numeric(time.1-sample, units = "mins")))

write.csv(tides_tidy, "tides.csv", row.names = F)

