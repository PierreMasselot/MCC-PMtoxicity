#####################################################################
#
#                 MCC-PMToxicity
#             Part 1: Preparing data
#
#####################################################################

#------------------------------------------
# Load MCC pollution data and keep cities with records
#------------------------------------------

#----- Load pollution data

# Load MCC data
mccpath <- "V:/VolumeQ/AGteam/MCCdata"
load(paste0(mccpath, "/air_pollution/MCC air pollution dataset/Processed", 
  "/MCCdata_Pollution_20230405.RData"))

# Select countries
cities <- subset(cities, country %in% country_select)
dlist <- dlist[cities$city]
countries <- subset(countries, country %in% country_select)

#----- Data cleaning

# Select observations
dlist <- lapply(dlist, subset, year >= yearstart)

# Loop
dlist <- lapply(dlist, mutate, 
  # Select outcome variable
  death = if(exists("all", mode = "numeric")) all else nonext,
  
  # Remove outliers
  death = if(exists("outlierm")) ifelse(outlierm == 1, NA, death) else death, 
  tmean = if(exists("outliert")) ifelse(outliert == 1, NA, tmean) else tmean,
  
  # Check for pm25
  pm25 = if(exists("pm25")) pm25 else NA,
  
  # Compute PM MA
  mapm = runMean(pm25, 0:maxlagp),
  
  # Convert O3 and NO2 to ppb and compute Ox
  no2 = if (exists('no2')) no2 / 1.88 else NA,
  o3 = if (exists('o3')) o3 / 1.96 else NA,
  ox = (no2 * oxw[1] + o3 * oxw[2]) / sum(oxw)
)

#----- Select data

# Select cities with enough data in common
haspm <- sapply(dlist, function(x) 
  sum(complete.cases(x[c("death", "mapm", "tmean")])) >= nday)
cities <- cities[haspm,]
dlist <- dlist[haspm]

#------------------------------------------
# Load city-level metadata
#------------------------------------------

#----- CAPI

# Load data
capidata <- read.csv("data/MCCpm25_locations_new.csv") |>
  na.omit()

#----- Ox data

# Compute average of oxidative pollutants
oxvars <- foreach(d = dlist, city = cities$city, .combine = rbind) %do% {
  summarise(d, city = city, o3_mcc = mean(o3, na.rm = T), 
    no2_mcc = mean(no2, na.rm = T), ox = mean(ox, na.rm = T))
} 


# #----- Compositions
# 
# # Load composition dataset
# flist <- list.files("data", pattern = "SPEC")
# years <- str_extract(flist, "[[:digit:]]{4}")
# pmcomp <- foreach(f = flist, y = years, .combine = rbind) %do% {
#   read.csv(sprintf("data/%s", f)) |> 
#     mutate(year = y)
# }
# 
# # Impute zero values
# pmcomp <- mutate(pmcomp, nodat = rowSums(pick(contains("PM25"))) == 0)
# impcomp <- dplyr::select(pmcomp, contains("PM25")) |> 
#   multRepl(label = 0, dl = rep(1e-5, 7), z.warning = 1)
# 
# # Create compositional data object
# comp_spec <- acomp(impcomp)
# 
# # Rename
# compnms <- names(comp_spec) |> strsplit("_") |> sapply("[", 2)
# names(comp_spec) <- compnms
# 
# # Average per city
# citycomp <- aggregate(comp_spec, subset(pmcomp, !nodat, city), mean)

#----- Load indicator data

# Load Urban Centre Database
load(paste0(mccpath, "/data/MCC_Indicators/MCC_indicators_UCD_20231110.RData"))

# Compute indicators of interest
ucd.mcc <- ucd.mcc |>
  mutate(GDPpc00 = GDP00_SM / P00, GDPpc15 = GDP15_SM / P15) |>
  dplyr::select(city, GDPpc00, GDPpc15, E_GR_AV00, E_GR_AV14, B00, B15)

#----- Temperature summary

# Compute mean and range temperature
tsum <- sapply(dlist, function(d) c(mean(d$tmean, na.rm = T), 
    diff(range(d$tmean, na.rm = T)))) |> 
  t() |> as.data.frame() |>
  rename(tmean = 1, trange = 2) |>
  rownames_to_column(var = "city")

#------------------------------------------
# Put everything together
#------------------------------------------

# Merge all metadata
cities <- Reduce(merge, list(cities, capidata, oxvars, ucd.mcc, tsum))

# Reorder
cities <- arrange(cities, country, city)

# Select series
dlist <- dlist[cities$city]
countries <- countries[countries$country %in% cities$country,]
