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
sensaiccs <- sapply(sensres, AICc)
sensbics <- sapply(sensres, BIC)

#----- Create final table

# Bind all results into table
sensrestab <- foreach(rer = sensrers, ci = sensrercis, p = senspvalues, 
  q = sensQs, i2 = sensi2s, lrt = senslrts, aic = sensaiccs, bic = sensbics, 
  lab = names(sensres), .combine = bind_rows) %do% 
{
  # RER string
  rerstr <- if (nrow(ci) > 0) sprintf("%.4f (%.4f - %.4f)", 
    rer, ci[,1], ci[,2]) else NULL
  
  # Create table
  out <- cbind(model = lab, var = names(rer), rer = rerstr, 
    lrt = lrt, aic = aic, bic = bic) |> as.data.frame()
  
}

# Column type
sensrestab <- mutate(sensrestab, across(all_of(c("lrt", "aic", "bic")), 
  as.numeric))


#----- Output into a word file
sensft <- flextable(sensrestab[, 
    c("model", "var", "rer", "lrt", "aic", "bic")]) |>
  # Merge identical rows
  merge_v(j = c("model", "lrt", "aic", "bic")) |>
  valign(j = c("model", "lrt", "aic", "bic"), valign = "top") |>
  # Put the best value in bold
  # bold(~ q == min(q), "q") |> 
  # bold(~ i2 == min(i2), "i2") |>
  bold(~ lrt == min(lrt, na.rm = T), "lrt") |>
  bold(~ aic == min(aic), "aic") |>
  bold(~ bic == min(bic), "bic") |>
  # Format columns
  colformat_double(j = c("aic", "bic"), digits = 2, big.mark = "") |>
  colformat_double(j = c("lrt"), digits = 4, na_str = "-") |>
  align(j = "var", align = "right")

# Change some labels
for (m in seq_along(modlabs)) sensft <- compose(sensft, 
  j = "model", i = ~ model == names(modlabs)[m], value = modlabs[[m]])
for (v in seq_along(varlabs)) sensft <- compose(sensft, 
  j = "var", i = ~ var == names(varlabs)[v], value = varlabs[[v]])


sensft <- sensft |>
  # Relabel
  set_header_labels(model = "Model", var = "", #q = "Cochran's Q", 
    lrt = "LRT P-value", aic = "AIC", bic = "BIC", rer = "RER (95% CI)") |>
  # Border
  fix_border_issues() |>
  # Resize
  autofit() |> fit_to_width(20, unit = "cm")
save_as_docx(sensft, path = "figures/SupTab6_sensitivityAnalysis.docx")

#----------------------------
# Comparison of RRs
#----------------------------

# List of featured countries
redcntr <- subset(cities, !is.na(coefadj), countryname, drop = T) |> unique()
senspal <- cntrpal[names(cntrpal) %in% redcntr]
sensshp <- shppal[names(shppal) %in% redcntr]

# Number of lower RRs
summarise(cities, lower = mean(coefadj < coef, na.rm = T), 
  meandiff = mean(exp(coef) - exp(coefadj), na.rm = T))
mutate(cities, diff = abs(exp(coef) - exp(coefadj))) |>
  arrange(desc(diff)) |>
  head()

# Plot
rrcompplot <- ggplot(cities) + theme_bw() + 
  geom_point(aes(x = exp(coef), y = exp(coefadj), 
    col = countryname, shape = countryname)) + 
  geom_hline(yintercept = 1) + 
  geom_vline(xintercept = 1) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values = senspal, breaks = names(senspal), 
    name = "Country", drop = T) + 
  scale_shape_manual(values = sensshp, breaks = names(senspal), 
    name = "Country", drop = T) + 
  labs(x = "RR (main)", y = "RR (adjusted)")
ggsave("figures/SupFig4_RRadj.png", rrcompplot, height = 5, width = 8)

#----------------------------
# Description of NO2 and O3 data
#----------------------------

#----- Summary by country

# Compute summaries of pollutants
senscntrsum <- subset(cities, convadj) |>
  summarise(ncities = n(), no2m = mean(NO2_ppbv), 
    no2q1 = quantile(NO2_ppbv, .25), no2q3 = quantile(NO2_ppbv, .75),
    o3m = mean(Ozone), o3q1 = quantile(Ozone, .25), 
    o3q3 = quantile(Ozone, .75),
    .by = countryname) |>
  arrange(countryname)

# Add total
senstotsum <- subset(cities, convadj) |>
  summarise(ncities = n(), no2m = mean(NO2_ppbv), 
    no2q1 = quantile(NO2_ppbv, .25), no2q3 = quantile(NO2_ppbv, .75),
    o3m = mean(Ozone), o3q1 = quantile(Ozone, .25), 
    o3q3 = quantile(Ozone, .75))
senscntrsum <- bind_rows(senscntrsum, senstotsum)

# Add name
senscntrsum <- merge(senscntrsum, 
  subset(countries, select = c(country, countryname)),
  all.x = T)
senscntrsum$countryname[is.na(senscntrsum$countryname)] <- "Total"

# Create pretty column
senscntrsum <- mutate(senscntrsum, 
  NO2 = sprintf("%.2f (%.2f - %.2f)", no2m, no2q1, no2q3),
  O3 = sprintf("%.2f (%.2f - %.2f)", o3m, o3q1, o3q3)) |>
  subset(select = c(countryname, ncities, NO2, O3))

#----- Create and export word table
senssumtab <- flextable(senscntrsum) |>
  set_header_labels(countryname = "Country", ncities = "Number of cities") |>
  compose(part = "header", j = "NO2", 
    value = as_paragraph("Average NO", as_sub("2"), " in ppbv (IQR)")) |>
  compose(part = "header", j = "O3", 
    value = as_paragraph("Average O", as_sub("3"), " in ppbv (IQR)")) |>
  bold(i = nrow(senscntrsum)) |>
  autofit() |> width(j = 2, width = .75)
save_as_docx(senssumtab, path = "figures/SupTab1_NO2O3_countryDesc.docx")