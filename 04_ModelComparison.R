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



