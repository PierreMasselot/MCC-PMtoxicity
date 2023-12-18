#####################################################################
#
#                 MCC-PMToxicity
#             Part 3: Second stage
#
#####################################################################

#----------------------------
# Control for city-characteristics: PCA
#----------------------------

# Run PCA
pcares <- prcomp(data.matrix(cities[, confvar]), center = T, scale = T)

# Extract components
pcs <- pcares$x[, 1:npc]

#----------------------------
# Main model: CAPI
#----------------------------

# Run model
capifit <- mixmeta(coef ~ PM2.5_toxicity + pcs, 
  S = v, random = ~ 1|country/city, data = cities, method = "ml", 
  subset = conv & !is.na(ox))
summary(capifit)

# Likelihood ratio tests
drop1(capifit, test = "Chisq")

#----------------------------
# Benchmark: Ox
#----------------------------

# Check Ox availability
sum(is.na(cities$ox))

# Fit model
oxfit <- mixmeta(coef ~ ox + pcs, 
  S = v, random = ~ 1|country/city, data = cities, 
  method = "ml", subset = conv & !is.na(ox))
summary(oxfit)

# Likelihood ratio tests
drop1(oxfit, test = "Chisq")

#----------------------------
# Sensitivity analysis 1: with adjusted RRs
#----------------------------

# Run model
adjfit <- mixmeta(coefadj ~ PM2.5_toxicity + pcs, 
  S = vadj, random = ~ 1|country/city, data = cities, method = "ml", 
  subset = convadj & !is.na(ox))
summary(adjfit)

# Likelihood ratio tests
drop1(adjfit, test = "Chisq")

#----------------------------
# Sensitivity analysis 1: with adjusted RRs by Ox
#----------------------------

# Run model
adjoxfit <- mixmeta(coefox ~ PM2.5_toxicity + pcs, 
  S = vox, random = ~ 1|country/city, data = cities, method = "ml", 
  subset = convox & !is.na(ox))
summary(adjoxfit)

# Likelihood ratio tests
drop1(adjoxfit, test = "Chisq")
