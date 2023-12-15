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
