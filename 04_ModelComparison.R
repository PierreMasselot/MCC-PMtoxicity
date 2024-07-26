#####################################################################
#
#                 MCC-PMToxicity
#             Part 4: Comparison with other effect modifiers
#
#####################################################################

# Object to include results
stage2res <- list()

# Add the main model
stage2res$PMCI <- pmcifit

#----------------------------
# Fit models for comparison
#----------------------------

#----- Null model

# Fit
stage2res$Null <- mixmeta(nullform, S = v, random = ranform, 
  data = cities, method = fitmethod, subset = conv)
summary(stage2res$Null)

# Extract main effect
predict(stage2res$Null, curvedata[1,, drop = F], ci = T) |> exp()

# Extract city-specific BLUPs
cities[cities$conv, "blup_null"] <- blup(stage2res$Null) |> exp()
summary(cities$blup_null)
subset(cities, rank(blup_null) <= 10, 
  c(cityname, countryname, PMCI, blup_null))
subset(cities, 
  rank(blup_null, na.last = F) >= (max(rank(blup_null, na.last = F)) - 10), 
  c(cityname, countryname, PMCI, blup_null))


#----- Model with gaseous pollutants

# List of variables
gasmod <- c("NO2_ppbv", "SO2", "Ozone", "HCHO", "CO", "NH3")

# Update formula
gasform <- update(nullform, 
  sprintf("~ . + %s", paste(gasmod, collapse = " + ")))

# Fit model
stage2res$gas <- mixmeta(gasform, S = v, random = ranform, 
  data = cities, method = fitmethod, subset = conv)
summary(stage2res$gas)

#----- Ox model

# Add Ox to formula
oxform <- update(nullform, ~ . + Ox)

# Fit model
stage2res$Ox <- mixmeta(oxform, S = v, random = ranform, 
  data = cities, method = fitmethod, subset = conv)
summary(stage2res$Ox)

#----- PM25 composition

# Transform components as log-ratios
alr_comp <- alr(cities[,comp_names])
p <- ncol(alr_comp)

# Update formula
compform <- update(nullform, ~ . + alr_comp)

# Fit model
stage2res$PM_Composition <- mixmeta(compform, S = v, random = ranform, 
  data = cities, method = fitmethod, subset = conv)
summary(stage2res$PM_Composition)
