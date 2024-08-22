#####################################################################
#
#                 MCC-PMToxicity
#            Second-stage functions
#
#####################################################################

#--------------------------
# PCA on confounders
#--------------------------

run_pca <- function(cities, confvar, npc = 1){
  
  # Run PCA
  pcares <- prcomp(data.matrix(cities[, confvar]), center = T, scale = T)
  
  # Extract components and return
  pcares$x[, 1:npc]
}

#--------------------------
# Meta-regression
#--------------------------

run_metareg <- function(cities, nullform, mainvar, ...){
  
  # Extract info on the main variable
  termlist <- attr(terms(mainvar), "term.label")
  
  #----- Fit model
  
  # Formula
  form <- reformulate(c(attr(terms(nullform), "term.label"), termlist),
    response = all.vars(nullform)[1])

  # Fit and display summary
  fit <- mixmeta(form, data = cities, ...)
  
  #----- Extract results
  
  # Extract coefficients and standard errors
  inds <- sapply(termlist, grep, names(coef(fit)), fixed = T) |> unlist()
  betas <- coef(fit)[inds]
  sebetas <-  sqrt(diag(vcov(fit)[inds, inds, drop = F]))
  
  # Compute RER for IQR increase with CI
  iqr <- apply(model.matrix(fit)[, inds, drop = F], 2, IQR) 
  rer <- exp(betas * iqr)
  rerci <- exp(iqr * (betas + 1.96 * t(c(-1, 1)) %x% sebetas))
  
  # Likelihood ratio test
  nullmod <- mixmeta(as.formula(deparse(nullform)), data = cities, ...)
  lrtres <- lrt.mixmeta(fit, nullmod)$pvalue
  
  # Information criteria
  ics <- c(AIC = AIC(fit), AICc = AICc(fit), BIC = BIC(fit))
  
  #----- Return
  
  # Rename variables
  varlist <- lapply(termlist, grep, colnames(model.matrix(fit)), 
    value = T, fixed = T) |> unlist()
  names(rer) <- rownames(rerci) <- varlist
  colnames(rerci) <- c("low", "high")
    
  # Return in list
  list(fit = fit, RER = cbind(rer, rerci), score = c(ics, LRT = lrtres))
}

#--------------------------
# Meta-regression for composition
#--------------------------

run_metareg_comp <- function(cities, nullform, mainvar, ...){
  
  # Extract variables and perform log-ratio transform
  termlist <- attr(terms(mainvar), "term.label")
  varlist <- all.vars(mainvar)
  
  #----- Fit model
  
  # Formula
  form <- reformulate(c(attr(terms(nullform), "term.label"), termlist),
    response = all.vars(nullform)[1])
  
  # Fit and display summary
  fit <- mixmeta(form, data = cities, ...)
  
  #----- Extract results
  
  # Extract coefficients and standard errors
  inds <- grep("alr", names(coef(fit)), fixed = T)
  betas <- coef(fit)[inds]
  sebetas <-  sqrt(diag(vcov(fit)[inds, inds, drop = F]))
  
  # Add for the last component
  betas <- c(betas, -sum(betas, na.rm = T))
  sebetas <- c(sebetas, sqrt(sum(vcov(fit)[inds, inds])))
  
  # Compute RER for IQR increase with CI
  iqr <- apply(cities[, varlist], 2, IQR) 
  rer <- exp(-2 * betas * log(iqr))
  rerci <- exp(-2 * log(iqr) * (betas + 1.96 * t(c(-1, 1)) %x% sebetas))
  
  # Likelihood ratio test
  nullmod <- mixmeta(as.formula(deparse(nullform)), data = cities, ...)
  lrtres <- lrt.mixmeta(fit, nullmod)$pvalue
  
  # Information criteria
  ics <- c(AIC = AIC(fit), AICc = AICc(fit), BIC = BIC(fit))
  
  #----- Return
  
  # Rename variables
  names(rer) <- rownames(rerci) <- varlist
  colnames(rerci) <- c("low", "high")
  
  # Return in list
  list(fit = fit, RER = cbind(rer, rerci), score = c(ics, LRT = lrtres))
}

#--------------------------
# Likelihood-ratio test
#--------------------------

lrt.mixmeta <- function(full, null) {
  
  # Extract log-likelihoods
  llnull <- logLik(null)
  llfull <- logLik(full)
  
  # Test statistic
  lrstat <- -2 * (llnull - llfull)
  attributes(lrstat) <- NULL
  
  # Compute p-value
  lrtdf <- attr(llfull, "df") - attr(llnull, "df")
  pval <- 1 - pchisq(lrstat, lrtdf)
  
  # Export test stat and p-value
  return(list(lrstat = lrstat, pvalue = pval))
}
