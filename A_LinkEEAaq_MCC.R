################################################################################
#
# Test linking EEAq data to MCC
#
################################################################################

library(EEAaq)
library(dplyr)
library(stringr)
library(utils)
library(foreach)
library(data.table)
library(lubridate)
library(tidyr)

#----- Parameters

# Pollutants
poll_list <- c("PM2.5", "O3") # c("PM10", "PM2.5", "O3", "NO2", "SO2")

# Years to extract pollutants
years <- c(2000, 2020)

# MCC dataset
mccdata <- "MCCdata_20230904"
mccpath <- "V:/VolumeQ/AGteam/MCCdata/data/MCC_all"

#---------------------------------
# Add EEAaq to MCC
#---------------------------------

#----- MCCPoll data

# Load
mccpollpath <- "V:/VolumeQ/AGteam/MCCdata"
load(paste0(mccpollpath, "/air_pollution/MCC air pollution dataset/Processed", 
  "/MCCdata_Pollution_20230405.RData"))

# Select MCC countries
country_poll <- c("aus8809", "can8615", "chi9615", "chl0414", "cyp0419", 
  "ecu1418", "est9718", "fnl9414", "fra0015", "ger9315", "grc0110", 
  "irn0215", "isr8520", "jap1115", "mex9814", "mlt9519", "nor6918", 
  "por8018", "rom9416", "sa9713", "spa9014", "sui9513", "swe9010", "twn9414",
  "uk9020", "usa7306")
cities_poll <- subset(cities, country %in% country_poll)
dlist_poll <- dlist[cities_poll$city]
countries_poll <- subset(countries, country %in% country_poll)

#----- MCC general dataset

# Load MCC and select countries
load(paste0(mccpath, "/", mccdata, ".RData"))

# Select countries
cntrsel <- c("bul0020", "cyp0419", "cze9420", "est9720", "fnl9414", "fra0017", 
  "ger9320", "grc0110", "ice0018", "irl8407", "ita0615", "mld0110", "mlt9519", 
  "net9516c", "nor6918", "por8018", "rom9416", "spa9014", "srb9521", "sui6918", 
  "swe9016", "uk9020")
cities_eu <- subset(cities, country %in% cntrsel)
dlist_eu <- dlist[cities_eu$city]
countries_eu <- subset(countries, country %in% cntrsel)

#----- Linking MCC and Eurostat

# Read lookup table
uraupath <- "V:/VolumeQ/AGteam/Eurostat/lookup"
mccurau <- read.table(paste0(uraupath, "/URAU_", mccdata, ".csv"), 
  sep = ";", quote = "\"", header = T)

# Select countries
mccurau <- mutate(mccurau, country = str_split_i(city, "\\.", 2)) |>
  subset(country %in% cntrsel)

# Remove those that are already in MCCpoll
mccurau <- subset(mccurau, !city %in% cities_poll$city)

#----- Download data

# Get list of stations operating in the time frame
eeastats <- EEAaq_get_stations() |> 
  subset(AirPollutant %in% poll_list &
    substr(OperationalActivityBegin, 7, 10) < years[2] & 
      (substr(OperationalActivityEnd, 7, 10) > years[1] | 
          is.na(OperationalActivityEnd)))

# Select LAUs that have at least one station for selected pollutants
laulist <- mutate(mccurau, LAU_ID = str_split_i(mccurau$LAU_CODE, "_", 2)) |>
  subset(select = c(LAU_ID, CNTR_CODE)) |>
  merge(unique(eeastats["LAU_ID"])) |>
  na.omit() |> unique()

# Download data for these LAUs
eupolldat <- foreach(lau = laulist$LAU_ID, cntr = laulist$CNTR_CODE,
  .final = rbindlist) %do% 
{
  print(lau)
  
  # Download (sometimes necessary to try several times)
  i <- 0
  repeat {
    dldat <- try(EEAaq_get_data(zone_name = lau, ID = T, NUTS_level = "LAU", 
      pollutant = poll_list, from = years[1], to = years[2], verbose = F))
    if(!inherits(dldat, "try-error") | i >= 2) break else i <- i + 1
  }
  
  # Check if it was downloaded properly
  if (inherits(dldat, "try-error")){
    return()
  } else {
    # Average hourly concentrations
    daydat <- mutate(dldat, date = date(with_tz(DatetimeBegin, tzone = "CET"))) |> 
      pivot_longer(cols = any_of(poll_list), names_to = "Pollutant") |>
      group_by(date, AirQualityStationEoICode, Pollutant) |>
      summarise(value = mean(value, na.rm = T), 
        complete = mean(is.na(value)) <= .25) |> 
      ungroup()
      
    # Clean and average stations
    mutate(daydat, value = ifelse(complete, value, NA)) |> 
      group_by(date, Pollutant) |>
      summarise(Concentration = mean(value, na.rm = T)) |>
      ungroup() |>
      mutate(LAU_CODE = paste(cntr, lau, sep = "_"), across())
  }
}

#----- Link to MCC cities

# Add city codes and aggregate by city
citypolldat <- merge(eupolldat, mccurau)[, 
  .(Concentration = mean(Concentration)), by = .(date, Pollutant, city)]

# Widen pollutants
citypolldat <- dcast(citypolldat, city + date ~ Pollutant)

# Add pollutants to cities
dlist_aug <- foreach(d = dlist_eu, cit = names(dlist_eu), 
  .final = function(x) {names(x) <- names(dlist_eu); x}) %do% 
{
  merge(d, subset(citypolldat, city == cit, -city), 
    by.x = "date", by.y = "Date", all.x = T)
}










# Keep only LAU-pollutant couples that exist
uni_laus <- unique(mccurau$LAU_CODE) |> na.omit()
dlgrid <- expand.grid(id = uni_laus, AirPollutant = poll_list) |>
  mutate(LAU_ID = str_split_i(id, "_", 2), 
    CNTR_CODE = str_split_i(id, "_", 1),
    id = NULL) |>
  merge(unique(eeastats[,c("LAU_ID", "AirPollutant")]))

# Download data by country and pollutant (to avoid timeout issues)
eupolldat <- foreach(x = split_laus, .final = rbindlist) %do% 
{
  # Download  
  dldat <- EEAaq_get_data(zone_name = x$LAU_ID, ID = T, NUTS_level = "LAU", 
    pollutant = as.character(unique(x$AirPollutant)), 
    from = years[1], to = years[2])
  
  # Aggregate to daily level
  EEAaq_time_aggregate(dldat, frequency = "daily", aggr_fun = "mean")
}






# Select stations with pollutants of interest
statlist <- subset(EEAaq::stations, AirPollutant %in% poll_list)

# Select stations with codes in lookup table
statlist <- mutate(statlist, LAU_CODE = paste(ISO, LAU_ID, sep = "_")) |>
  subset(LAU_CODE %in% mccurau$LAU_CODE)

# Download
eupolldat <- lapply(poll_list, function(p) 
  EEAaq_get_data(zone_name = na.omit(statlist$LAU_ID), ID = T, 
    NUTS_level = "LAU", pollutant = p, 
    from = years[1], to = years[2]))

#----- Aggregate for MCC cities

# Aggregate to the daily level
eupolldat <- EEAaq_time_aggregate(eupolldat, frequency = "daily", 
  aggr_fun = "mean")
if (length(poll_list) > 1) {
  eupolldat <- eupolldat[[1]]
} else {
  eupolldat <- rename_with(eupolldat, 
    ~ paste0(poll_list, "_mean"), ends_with("mean"))
}
names(eupolldat) <- gsub("_mean", "", names(eupolldat))

# link station to MCC cities
eupolldat <- as.data.table(eupolldat)
eupolldat <- merge(eupolldat, 
    unique(statlist[,c("AirQualityStationEoICode", "LAU_CODE")]),
    by = "AirQualityStationEoICode") |>
  merge(na.omit(mccurau[,c("LAU_CODE", "city")]), by = "LAU_CODE", 
    allow.cartesian = T)

# Aggregate pollutants for each city
eupolldat <- eupolldat[, lapply(.SD, mean), by = .(city, Date), 
  .SDcols = patterns(paste(poll_list, collapse = "|"))]

#----- Add to MCC dlist

# Load MCC and select countries
load(paste0(mccpath, "/", mccdata, ".RData"))
cities_eu <- subset(cities, country %in% cntrsel)
dlist_eu <- dlist[cities_eu$city]
countries_eu <- subset(countries, country %in% cntrsel)

# Add pollutants to cities
dlist_eu <- foreach(d = dlist_eu, cit = names(dlist_eu), 
  .final = function(x) {names(x) <- names(dlist_eu); x}) %do% 
{
  merge(d, subset(eupolldat, city == cit, -city), by.x = "date", by.y = "Date",
    all.x = T)
}


#---------------------------------
# Compare with MCCpoll
#---------------------------------

#----- Load MCC poll

# Load
mccpollpath <- "V:/VolumeQ/AGteam/MCCdata"
load(paste0(mccpollpath, "/air_pollution/MCC air pollution dataset/Processed", 
  "/MCCdata_Pollution_20230405.RData"))

# Select MCC countries
country_poll <- c("aus8809", "can8615", "chi9615", "chl0414", "cyp0419", 
  "ecu1418", "est9718", "fnl9414", "fra0015", "ger9315", "grc0110", 
  "irn0215", "isr8520", "jap1115", "mex9814", "mlt9519", "nor6918", 
  "por8018", "rom9416", "sa9713", "spa9014", "sui9513", "swe9010", "twn9414",
  "uk9020", "usa7306")
cities_poll <- subset(cities, country %in% country_poll)
dlist_poll <- dlist[cities_poll$city]
countries_poll <- subset(countries, country %in% country_poll)

#----- Compare

# List of EU cities with at least one year of pollution
haspoll <- lapply(dlist_eu, function(d) summarise(d, across(all_of(poll_list), 
  ~ sum(!is.na(.x)) >= 365))) |> rbindlist()
pollcities <- lapply(haspoll, function(x) names(dlist_eu)[x])

# Check which are not in MCC poll
diffcity <- lapply(pollcities, setdiff, cities_poll$city)

# Same but with names
diffcityname <- lapply(pollcities, function(x) 
  subset(cities_eu, city %in% x, c("cityname", "countryname")) |>
    setdiff(cities_poll[,c("cityname", "countryname")]))
sapply(diffcityname, NROW)