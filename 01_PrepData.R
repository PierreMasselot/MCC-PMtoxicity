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
  "/MCCdata_Pollution_20240313.RData"))

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
  
  # Check for daily gases
  no2 = if (exists('no2')) no2 else NA,
  o3 = if (exists('o3')) o3 else NA
)

#----- Select data

# Select cities with enough data in common
haspm <- sapply(dlist, function(x) 
  sum(complete.cases(x[c("death", "mapm", "tmean")])) >= nday)
cities <- cities[haspm,]
dlist <- dlist[haspm]

#----- Separate USA in regions

# Load data for regions
path <- "V:/VolumeQ/AGteam/MCCdata/original/USA211cities/Final_Dataset"
load(paste(path, "usa7306.RData",sep="/"))

# Match regions to cities
usaregion <- data.frame(Region=1:9, regionname=paste("USA", c("Central", 
  "NECentral", "NWCentral", "NorthEast", "NorthWest", "South", "SouthEast",
  "SouthWest", "West"), sep="-"))
temp <- merge(citiesusa7306[c("city","Region")], usaregion)

# Replace in meta dataset
ind <- match(temp$city, cities$city)
cities[na.omit(ind), "country"] <- cities[na.omit(ind), "countryname"] <- 
  temp$regionname[!is.na(ind)]


#------------------------------------------
# Load city-level metadata
#------------------------------------------

#----- CAPI data

# Load data
capiread <- read.csv("data/city_pollutants_v2.csv")

# Compute Ox
polldata <- mutate(capiread, Ox = (NO2_ppbv * 1.07 + Ozone * 2.075) / 3.145)

# Remove country information
polldata <- subset(polldata, !is.na(PMCI), select = -c(country, countryname))

#----- PM composition data

# Find files
pathcomp <- paste0(mccpath, "/air_pollution/Composition_data/MCC_PM_SPEC_V2")
compfiles <- list.files(pathcomp, pattern = "SPEC_10km_buffer")
compyears <- str_split_i(compfiles, pattern = "[_\\.]", 6)

# Load all files
pmcomp <- Map(function(f, y) read.csv(paste0(pathcomp, "/", f)) |> 
    mutate(year = y),
  compfiles, compyears)
pmcomp <- do.call(rbind, pmcomp)

# Variables
comp_inds <- grep("PM25", colnames(pmcomp))
comp_names <- str_split_i(colnames(pmcomp)[comp_inds], "_", 2)

# Loop on cities to compute average composition
compmeans <- split(pmcomp, ~ city) |> sapply(function(d){
  
  # Zeros imputation
  allzeros <- apply(d[,comp_inds], 1, function(x) all(x == 0))
  imp_comp <- multRepl(d[!allzeros, comp_inds], label = 0, dl = rep(1e-5, 7),
    z.warning = 1.1)
  
  # Compositional object and mean
  acomp(imp_comp) |> mean()
})

# Tidy
compmeans <- t(compmeans) |> as.data.frame()
colnames(compmeans) <- comp_names
compmeans$city <- rownames(compmeans)

#----- Load indicator data

# Load Urban Centre Database
load(paste0(mccpath, "/data/MCC_Indicators/MCC_indicators_UCD_20240313.RData"))

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
cities <- Reduce(merge, list(cities, polldata, compmeans, ucd.mcc, tsum))

# Reorder
cities <- arrange(cities, country, city)

# Select series
dlist <- dlist[cities$city]
countries <- countries[countries$country %in% cities$country,]
