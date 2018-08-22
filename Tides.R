library(dplyr)
library(rvest)
library(readr)

# station IDs

## Whidbey
# 644 <- CrescentHarbor (MarinersCove, CrescentHarbor, ForbesPoint, MaylorPoint)
# 620 <- Coupeville (MonroeLanding, RollingHills, Coupeville)
# 1078 <- Greenbank (HarringtonN, HarringtonS, PrattsBluff)
# 2395 <- SandyPoint (Langley)
# 1026 <- Glendale (PossessionPoint)
# 361 <- BushPoint (Bush Point Dock, MalmoBluff, ShoreMeadows, MutinySands, LimpetLane, DoubleBluff
#                   LagoonN, LagoonS, Keystone, Ledgewood, Hancock, FortCasey)
# 2686 <- SunsetBeach (Cliffside, Swantown, HastieLake)

## South Sound
# 155 <- Totten Legacy, Totten Elizan C, Totten Elizan A Tower, Totten Elizan A North, 
# 2306 <- Youngs Cove A, Youngs Cove B, Youngs Cove C, Youngs Cove D, Youngs Cove E, Youngs Cove F, Flapjack
# 1909 <- Priest Point, Hearthfire
# 344 <- Gull Harbor, Burfoot Park, Burfoot Extended
# 733 <- Zangle Cove A, Big Fish Trap, Briscoe Pt
# 715 <- Amsterdam Bay, Mill Bight A, Mill Bight B, Mill Bight Extended, Walnut Rd A, Walnut Rd B, Walnut Rd C, Beachcrest,
#           Butterball Cove South A, Butterball North, Andys Marine Park South, Andys Marine Park West, Lyle Point
# 3048 <- Sandy Pt, Higgins Cove
# 2643 <- Ketron Ferry, Ketron SE, Ketron SW
# 2906 <- Jarrell Cove State Park
# 773 <- Cole Pt
# 1147 <- Zittels Marina


col_id <- data.frame(site = 
        c("Crescent Harbor", "Forbes Point", "Maylor Point",
          "Coupeville Wharf", "Monroe Landing", "Rolling Hills",
          "Harrington North", "Harrington South", "Pratts Bluff",
          "Langley Marina",
          "Possession Point",
          "Bush Point Dock", "Malmo Bluff", "Shore Meadows", "Mutiny Sands", "Limpet Lane", 
                "Double Bluff North", "Double Bluff South", "Lagoon North #0", "Lagoon North #1", 
                "Lagoon North #2", "Lagoon North #3", "Lagoon South", 
                "Keystone", "Ledgewood", "Hancock Lake", "Fort Casey", "Libby North",
          "Cliffside", "Swantown", "Hastie Lake", "Hastie Lake North"), 
                             grp = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, rep(6, 17),
                                     7, 7, 7, 7))

station_id <- data.frame(station = c(644, 620, 1078, 2395, 1026, 361, 2686), grp = 1:7) %>%
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
  transform(sample = as.POSIXlt(paste(date, 
              format(strptime("8:00", format = "%H:%M"), "%H:%M:%S"), sep = " "))) %>%
  transform(from_low = 
              ifelse(type == "L", as.numeric(time-sample, units = "mins"), 
                     as.numeric(time.1-sample, units = "mins")),
            from_high = 
              ifelse(type == "H", as.numeric(time-sample, units = "mins"), 
                     as.numeric(time.1-sample, units = "mins"))) %>%
  merge(station_id %>% select(-grp), by = 'station')

write.csv(tides_tidy, "tides_tidy.csv", row.names = F)

