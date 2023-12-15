#####################################################################
#
#                 MCC-PMToxicity
#             Test control for city-level variables
#
#####################################################################

# List of confounders
confvar <- c("GDPpc00", "GDPpc15", "E_GR_AV00", "E_GR_AV14", "B00", 
  "B15", "tmean", "trange", "PM2.5_Linear")
nconf <- length(confvar)

#-------------------------
# PCA
#-------------------------

# Run PCA
pcares <- prcomp(data.matrix(cities[, confvar]), center = T, scale = T)
summary(pcares)

# Extract components
pcavar <- pcares$x

# Screeplot
var_explained <- pcares$sdev^2 / sum(pcares$sdev^2)
plot(var_explained, type = "b", xlab = "Components", 
  ylab = "Variance proportion", ylim = c(0, 1), pch = 16)
points(cumsum(var_explained), type = "b", pch = 16, col = 2)


#-------------------------
# PLS
#-------------------------

# Run PLS
plsform <- sprintf("PM2.5_toxicity ~ %s", paste(confvar, collapse = " + "))
plsres <- plsr(as.formula(plsform), data = cities, scale = T)

# Extract scores for all cities
plsvar <- predict(plsres, newdata = cities, type = "scores")
colnames(plsvar) <- sprintf("pls%i", seq_len(ncol(plsvar)))

#-------------------------
# Propensity function
#-------------------------

# Run linear model
lmform <- sprintf("PM2.5_toxicity ~ %s", paste(confvar, collapse = " + "))
lmres <- lm(as.formula(lmform), data = cities)

# Extract prediction (PF)
pfun <- predict(lmres)


#-------------------------
# Compare balance
#-------------------------

#----- No control

basebal <- foreach(v = confvar, .combine = c) %do% 
  {
    # Regression on covariate
    form <- sprintf("%s ~ PM2.5_toxicity", v)
    regres <- lm(form, cities)
    
    # Extract t-value
    summary(regres)$coefficients["PM2.5_toxicity", "t value"]
  }

#----- PCA

# Loop on components
pcabal <- foreach(npc = seq_along(confvar), .combine = cbind) %:% 
  foreach(v = confvar, .combine = c) %do% 
{
    
  # Regression on covariate
  form <- sprintf("%s ~ PM2.5_toxicity + pc", v)
  pc <- pcavar[, 1:npc]
  regres <- lm(form, cities)
  
  # Extract t-value
  summary(regres)$coefficients["PM2.5_toxicity", "t value"]
}

# Plot T distribution
boxplot(abs(pcabal), xlab = "Components", ylab = "T statistics", 
  names = seq_along(confvar))


#----- PLS

# Loop on components
plsbal <- foreach(npc = seq_along(confvar), .combine = cbind) %:% 
  foreach(v = confvar, .combine = c) %do% 
  {
    
    # Regression on covariate
    form <- sprintf("%s ~ PM2.5_toxicity + pls", v)
    pls <- plsvar[, 1:npc]
    regres <- lm(form, cities)
    
    # Extract t-value
    summary(regres)$coefficients["PM2.5_toxicity", "t value"]
  }

# Plot T distribution
boxplot(abs(plsbal), xlab = "Components", ylab = "T statistics", 
  names = seq_along(confvar))

#----- PF

pfbal <- foreach(v = confvar, .combine = c) %do% 
{
  # Regression on covariate
  form <- sprintf("%s ~ PM2.5_toxicity + pfun", v)
  regres <- lm(form, cities)
  
  # Extract t-value
  summary(regres)$coefficients["PM2.5_toxicity", "t value"]
}

#----- Compare everything

# Put everything together
allbal <- cbind(basebal, pcabal, plsbal, pfbal)
balnames <- c("None", sprintf("PCA%i", 1:nconf), sprintf("PLS%i", 1:nconf), 
  "PFlm")

# Compute means and cis
means <- apply(abs(allbal), 2, mean)
lows <- apply(abs(allbal), 2, quantile, .025)
highs <- apply(abs(allbal), 2, quantile, .975)

# Plot
plot(0, type = "n", xlim = c(1, ncol(allbal)), ylim = range(c(lows, highs)),
  xlab = "", ylab = "T-statistic distribution", xaxt = "n")
arrows(x0 = 1:ncol(allbal), y0 = lows, y1 = highs, len = .05, angle = 90, 
  code = 3)
points(means, pch = 16)
axis(1, at = 1:ncol(allbal), labels = balnames, las = 2)


#-------------------------
# Compare obtained effect
#-------------------------

#----- No control

# Fit mixmeta
mixbase <- mixmeta(coef ~ PM2.5_toxicity, 
  S = v, random = ~ 1|country/city, data = cities, method = "ml", 
  subset = conv)

# Extract effect
basecoef <- summary(mixbase)$coefficients["PM2.5_toxicity", 
  c("Estimate", "95%ci.lb", "95%ci.ub")]

#----- PCA

# Loop on components
pcacoef <- foreach(npc = seq_along(confvar), .combine = cbind) %do% 
  {
    
    # Regression on covariate
    pc <- pcavar[, 1:npc]
    mixres <- mixmeta(coef ~ PM2.5_toxicity + pc, 
      S = v, random = ~ 1|country/city, data = cities, method = "ml", 
      subset = conv)
    
    # Extract t-value
    summary(mixres)$coefficients["PM2.5_toxicity", 
      c("Estimate", "95%ci.lb", "95%ci.ub")]
  }

#----- PLS

# Loop on components
plscoef <- foreach(npc = seq_along(confvar), .combine = cbind) %do% 
  {
    
    # Regression on covariate
    pls <- plsvar[, 1:npc]
    mixres <- mixmeta(coef ~ PM2.5_toxicity + pls, 
      S = v, random = ~ 1|country/city, data = cities, method = "ml", 
      subset = conv)
    
    # Extract t-value
    summary(mixres)$coefficients["PM2.5_toxicity", 
      c("Estimate", "95%ci.lb", "95%ci.ub")]
  }

#----- PF

# Fit mixmeta
mixpf <- mixmeta(coef ~ PM2.5_toxicity + pfun, 
  S = v, random = ~ 1|country/city, data = cities, method = "ml", 
  subset = conv)

# Extract effect
pfcoef <- summary(mixpf)$coefficients["PM2.5_toxicity", 
  c("Estimate", "95%ci.lb", "95%ci.ub")]

#----- Compare everything

# Put everything together
allcoef <- cbind(basecoef, pcacoef, plscoef, pfcoef)
coefnames <- c("None", sprintf("PCA%i", 1:nconf), sprintf("PLS%i", 1:nconf), 
  "PFlm")

# Plot
plot(0, type = "n", xlim = c(1, ncol(allcoef)), ylim = range(c(allcoef, 0)),
  xlab = "", ylab = "PM2.5 Toxicity coefficient", xaxt = "n")
arrows(x0 = 1:ncol(allcoef), y0 = allcoef[3,], y1 = allcoef[2,], 
  len = .05, angle = 90, code = 3)
points(allcoef[1,], pch = 16)
axis(1, at = 1:ncol(allcoef), labels = coefnames, las = 2)
abline(h = 0)
