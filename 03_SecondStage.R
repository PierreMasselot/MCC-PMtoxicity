#####################################################################
#
#                 MCC-PMToxicity
#             Part 3: Second stage
#
#####################################################################

#----------------------------
# Fit meta-regression model with Toxicity Index
#----------------------------

#----- PCA on confounders

# Run PCA
pcares <- prcomp(data.matrix(cities[, confvar]), center = T, scale = T)

# Extract components
pcs <- pcares$x[, 1:npc]

#----- Meta-regression model

# Formula
pmciform <- update(nullform, ~ . + PMCI)

# Fit and display summary
pmcifit <- mixmeta(pmciform, S = v, random = ranform, 
  data = cities, method = fitmethod, subset = conv)
summary(pmcifit)

#----------------------------
# Predict RR
#----------------------------

#----- Extract BLUP and residuals

# BLUPs
cities[cities$conv, "blup"] <- exp(blup(pmcifit))

# Residuals
cities[cities$conv, "residuals"] <- residuals(pmcifit)

#----- Predict RR for a grid of PMCI

# Create grid of PMCI
pmcigrid <- with(cities, seq(min(PMCI), max(PMCI), length.out = 100))
pcmat <- matrix(0, nrow = 100, ncol = 2, dimnames = list(NULL, colnames(pcs)))
newdata <- data.frame(PMCI = pmcigrid, pcs = I(pcmat))

# Predict
rrpred <- predict(pmcifit, newdata, ci = T) |> 
  exp() |>
  as.data.frame()
newdata <- cbind(newdata, rrpred)

#----- Predict RR at IQR

# Quantile to predict
predper <- c(.25, .75)

# Newdata
perdata <- data.frame(PMCI = quantile(cities$PMCI, predper), 
  pcs = I(pcmat[seq_along(predper),]))

# Predict
predict(pmcifit, perdata, ci = T) |> exp()

