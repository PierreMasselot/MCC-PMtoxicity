#####################################################################
#
#                 MCC-PMToxicity
#             Part 2: First stage, city level
#
#####################################################################

#-------------------------------------
# Loop on cities
#-------------------------------------

# Prepare parallelisation
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Loop
stage1res <- foreach(dat = dlist, city = names(dlist), 
  .packages = c("dlnm", "splines"), .combine = rbind) %dopar% 
{

  # Construct crossbasis for temperature confounding
  cbt <- crossbasis(dat$tmean, lag = 3, 
    arglag = list(fun = "strata"),
    argvar = list(fun = "bs", 
      knots = quantile(dat$tmean, c(.1, .75, .9), na.rm = T))
  )
  
  # Construct crossbasis for PM2.5
  cbp <- crossbasis(dat$pm25, lag = maxlagp, argvar = argvarp, arglag = arglagp) 
  
  # Estimate the model
  model <- glm(death ~ cbp + cbt + dow + 
      ns(date, df = timedf * length(unique(year))), 
    data = dat, family = quasipoisson)
  
  # Extract logRR for an increase of 10ug
  redall <- if (model$converged){
    crosspred(cbp, model, cen = cen, at = 10)
  } else {
    list(allfit = NA, allse = NA)
  }
  
  # Return
  data.frame(city = city, coef = redall$allfit, v = redall$allse^2, 
    conv = model$converged)
}
  

# Close parallel
stopCluster(cl)

# Add to cities object
cities <- merge(cities, stage1res)

# Save intermediate results
save.image("temp/fs_res.RData")