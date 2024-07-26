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
pmciform <- update(nullform, ~ . + log(I(PMCI + 1)))

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
cities[cities$conv, "mar_resid"] <- residuals(pmcifit)
cities[cities$conv, "cntr_resid"] <- blup(pmcifit, type = "residual", level = 1)

#----- Predict RR for a grid of PMCI

# Create grid of PMCI
pmcigrid <- with(cities, seq(min(PMCI), max(PMCI), length.out = 100))
pcmat <- matrix(0, nrow = 100, ncol = 2, dimnames = list(NULL, colnames(pcs)))
curvedata <- data.frame(PMCI = pmcigrid, pcs = I(pcmat))

# Predict
rrpred <- predict(pmcifit, curvedata, ci = T) |> 
  exp() |>
  as.data.frame()
curvedata <- cbind(curvedata, rrpred)

#----- Predict RR at IQR

# Quantile to predict
predper <- c(.25, .75)

# Data for prediction
perdata <- data.frame(PMCI = quantile(cities$PMCI, predper), 
  pcs = I(pcmat[seq_along(predper),]))

# Predict
predict(pmcifit, perdata, ci = T) |> exp() |> as.data.frame() |>
  reframe(RRres = sprintf("%.4f (95%%CI: %.4f - %.4f)", fit, ci.lb, ci.ub))
