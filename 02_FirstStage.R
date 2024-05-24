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
stage1res <- foreach(dat = dlist, city = iter(cities, by = "row"), 
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
  
  #----- Main model
  
  # Estimate the model
  model <- glm(death ~ cbp + cbt + dow + 
      ns(date, df = timedf * length(unique(year))), 
    data = dat, family = quasipoisson)
  
  # Extract logRR for an increase of 10ug
  redall <- if (model$converged){
    crosspred(cbp, model, cen = 0, at = 10)
  } else {
    list(allfit = NA, allse = NA)
  }
  
  #----- Adjustment for O3 and NO2
  
  if (!(is.null(dat$o3) | is.null(dat$no2) | 
      all(is.na(dat$o3)) | all(is.na(dat$no2))))
  {
    
    # Crossbasis for O3
    cbo3 <- crossbasis(dat$o3, lag = 1, argvar = list(fun = "lin"), 
      arglag = list(fun = "strata")) 
    
    # Crossbasis for NO2
    cbno2 <- crossbasis(dat$no2, lag = 1, argvar = list(fun = "lin"), 
      arglag = list(fun = "strata")) 
    
    # Estimate model
    modeladj <- glm(death ~ cbp + cbo3 + cbno2 + cbt + dow + 
        ns(date, df = timedf * length(unique(year))), 
      data = dat, family = quasipoisson)
    
    # Extrac logRR
    redalladj <- crosspred(cbp, modeladj, cen = 0, at = 10)
    convadj <- modeladj$converged
      
  } else {
    redalladj <- list(allfit = NA, allse = NA)
    convadj <- F
  }
  
  #----- Return results
  data.frame(city = city$city,
    coef = redall$allfit, v = redall$allse^2, conv = model$converged, 
    coefadj = redalladj$allfit, vadj = redalladj$allse^2, convadj = convadj,
    deaths = sum(model.frame(model)$death),
    periodmin = dat[-attr(model.frame(model), "na.action"), "year"] |> min(), 
    periodmax = dat[-attr(model.frame(model), "na.action"), "year"] |> max())
}

# Close parallel
stopCluster(cl)

# Add to cities object
cities <- merge(cities, stage1res, sort = F)
