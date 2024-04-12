#####################################################################
#
#                 MCC HetPoll
#               Part 6: Sensitivity analysis
#
#####################################################################

# Object to include results
sensres <- list()

# Change the null formula
sensform <- update(nullform, coefadj ~ .)

#----------------------------
# Fit models for comparison
#----------------------------

# Main model
senspmciform <- update(sensform, ~ . + log(I(PMCI + 1)))
sensres$PMCI <- mixmeta(senspmciform, S = vadj, random = ranform, 
  data = cities, method = fitmethod, subset = convadj)

# Null model
sensres$Null <- mixmeta(sensform, S = vadj, random = ranform, 
  data = cities, method = fitmethod, subset = convadj)

# Ox model
sensoxform <- update(sensform, ~ . + Ox)
sensres$Ox <- mixmeta(sensoxform, S = vadj, random = ranform, 
  data = cities, method = fitmethod, subset = convadj)

# PM25 composition
senscompform <- update(sensform, ~ . + alr_comp)
sensres$PM_Composition <- mixmeta(senscompform, S = vadj, random = ranform, 
  data = cities, method = fitmethod, subset = convadj)

#----------------------------
# Table of results
#----------------------------

#----- Extract effect modification results

# Compute IQRs
sensiqrs <- Map(function(x, v) apply(model.matrix(x)[, v, drop = F], 2, IQR),
  sensres, vars)

# Extract coefficients
sensmetacoefs <- Map(function(x, v) coef(x)[v], sensres, vars)

# Extract standard errors
sensses <- Map(function(x, v) sqrt(diag(vcov(x)[v, v, drop = F])), 
  sensres, vars)

# Add last component for PM composition
sensiqrs$PM_Composition <- -log(c(alrInv(sensiqrs$PM_Composition, 
  orig = alr_comp)))
sensmetacoefs$PM_Composition <- c(sensmetacoefs$PM_Composition, 
  -sum(sensmetacoefs$PM_Composition))
sensses$PM_Composition <- c(sensses$PM_Composition, 
  sqrt(sum(vcov(sensres$PM_Composition)[vars$PM_Composition, 
    vars$PM_Composition])))
names(sensiqrs$PM_Composition) <- names(sensmetacoefs$PM_Composition) <- 
  names(sensses$PM_Composition) <- comp_names

# Compute RERs for IQR increase
sensrers <- Map(function(co, iqr) exp(co * iqr),
  sensmetacoefs, sensiqrs)

# Compute confidence intervals
sensrercis <- Map(function(co, s, iqr) 
  exp(cbind(co - 1.96 * s, co + 1.96 * s) * iqr),
  sensmetacoefs, sensses, sensiqrs)

# Extract p-value
senspvalues <- Map(function(co, s) 2 * (1 - pnorm(abs(co / s))),
  sensmetacoefs, ses)

#----- Extract model criteria

# Extract Cochran's Q and I2
sensQs <- sapply(sensres, function(x) summary(x)$qstat$Q)
sensi2s <- sapply(sensres, function(x) summary(x)$i2stat)

# Extract LRT p-value
senslrts <- sapply(sensres, function(x) lrt.mixmeta(x, sensres$Null)$pvalue)

# Extract AICs
sensaics <- sapply(sensres, AIC)

#----- Create final table

# Labels
modlabs <- c("Main", "Null", "Ox", "PM Composition")

# Bind all results into table
sensrestab <- foreach(rer = sensrers, ci = sensrercis, p = senspvalues, 
  q = sensQs, i2 = sensi2s, lrt = senslrts, aic = sensaics, lab = modlabs, 
  .combine = bind_rows) %do% 
{
  # RER string
  rerstr <- if (nrow(ci) > 0) sprintf("%.4f (%.4f - %.4f)", 
    rer, ci[,1], ci[,2]) else NULL
  
  # Create table
  out <- cbind(model = lab, var = names(rer), rer = rerstr, 
    q = q, i2 = i2, lrt = lrt, aic = aic) |> as.data.frame()
  
}

# Column type
sensrestab <- mutate(sensrestab, across(all_of(c("q", "i2", "lrt", "aic")), 
  as.numeric))


#----- Output into a word file
sensft <- flextable(sensrestab[, 
  c("model", "var", "rer", "q", "i2", "lrt", "aic")]) |>
  # Merge identical rows
  merge_v(j = c("model", "q", "i2", "lrt", "aic")) |>
  valign(j = c("model", "q", "i2", "lrt", "aic"), valign = "top") |>
  # Put the best value in bold
  bold(~ q == min(q), "q") |> 
  bold(~ i2 == min(i2), "i2") |>
  bold(~ lrt == min(lrt, na.rm = T), "lrt") |>
  bold(~ aic == min(aic), "aic") |>
  # Format columns
  colformat_double(j = c("q", "i2", "aic"), digits = 0, big.mark = "") |>
  colformat_double(j = c("lrt"), digits = 4, na_str = "-") |>
  align(j = "var", align = "right") |>
  # Relabel header
  set_header_labels(model = "Model", var = "", q = "Cochran's Q", 
    lrt = "LRT P-value", aic = "AIC", rer = "RER (95% CI)") |>
  compose(part = "header", j = "i2", value = as_paragraph("I", as_sup("2"))) |>
  # Border
  fix_border_issues() |>
  # Resize
  autofit()
save_as_docx(sensft, path = "figures/SupTab1_sensitivityAnalysis.docx")

