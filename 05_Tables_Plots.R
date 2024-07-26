#####################################################################
#
#                 MCC HetPoll
#               Part 5: Tables and Plots
#
#####################################################################

#-----------------------
# Table 1: Description of cities
#-----------------------

#----- Summary by country

# Compute summaries of pollutants
cntrsum <- subset(cities, conv) |>
  summarise(ncities = n(), PM25m = mean(PM25), 
    PM25q1 = quantile(PM25, .25), PM25q3 = quantile(PM25, .75),
    pmcim = mean(PMCI), pmciq1 = quantile(PMCI, .25), 
    pmciq3 = quantile(PMCI, .75),
    deaths = sum(deaths), 
    periodmin = min(periodmin), periodmax = max(periodmax),
  .by = countryname) |>
  arrange(countryname)

# Add total
totsum <- subset(cities, conv) |>
  summarise(ncities = n(), PM25m = mean(PM25), 
    PM25q1 = quantile(PM25, .25), PM25q3 = quantile(PM25, .75),
    pmcim = mean(PMCI),
    pmciq1 = quantile(PMCI, .25), 
    pmciq3 = quantile(PMCI, .75),
    deaths = sum(deaths), 
    periodmin = min(periodmin), periodmax = max(periodmax))
cntrsum <- bind_rows(cntrsum, totsum)

# Add name
cntrsum <- merge(cntrsum, subset(countries, select = c(country, countryname)),
  all.x = T)
cntrsum$countryname[is.na(cntrsum$countryname)] <- "Total"

# Create pretty column
cntrsum <- mutate(cntrsum, 
    PM25 = sprintf("%.2f (%.2f - %.2f)", PM25m, PM25q1, PM25q3),
    PMCI = sprintf("%.2f (%.2f - %.2f)", pmcim, pmciq1, pmciq3),
    period = sprintf("%i - %i", periodmin, periodmax)) |>
  subset(select = c(countryname, ncities, period, deaths, PM25, PMCI))

#----- Create and export word table
sumtab <- flextable(cntrsum) |>
  set_header_labels(countryname = "Country", ncities = "Number of cities",
    period = "Period", deaths = "Total mortality", 
    PMCI = "Average PMCI (IQR)") |>
  compose(part = "header", j = 5, 
    value = as_paragraph("Average PM", as_sub("2.5"), " in µg/m", as_sup("3"), " (IQR)")) |>
  bold(i = nrow(cntrsum)) |>
  autofit() |> #width(j = 2, width = .75) |> width(j = 4, width = 1) |>
  fit_to_width(20, unit = "cm")

# Save
save_as_docx(sumtab, path = "figures/Tab1_countryDesc.docx")

#----- Additional results

# Percentage locations with average below 10ug/m3
subset(cities, conv) |> summarise(mean(PM25 < 10) * 100)

# Percentage locations with positive PMCI
subset(cities, conv) |> summarise(mean(PMCI > 0) * 100)

#-----------------------
# Table 2: Model result and comparison
#-----------------------

#----- Extract effect modification results

# List of variables to extract
nullvars <- model.matrix(stage2res$Null) |> colnames()
vars <- lapply(stage2res, function(x) setdiff(colnames(model.matrix(x)),
  nullvars))

# Compute IQRs
iqrs <- Map(function(x, v) apply(model.matrix(x)[, v, drop = F], 2, IQR),
  stage2res, vars)

# Extract coefficients
metacoefs <- Map(function(x, v) coef(x)[v], stage2res, vars)

# Extract standard errors
ses <- Map(function(x, v) sqrt(diag(vcov(x)[v, v, drop = F])), stage2res, vars)

# Add last component for PM composition
iqrs$PM_Composition <- -log(c(alrInv(iqrs$PM_Composition, orig = alr_comp)))
metacoefs$PM_Composition <- c(metacoefs$PM_Composition, 
  -sum(metacoefs$PM_Composition))
ses$PM_Composition <- c(ses$PM_Composition, 
  sqrt(sum(vcov(stage2res$PM_Composition)[vars$PM_Composition, 
    vars$PM_Composition])))
names(iqrs$PM_Composition) <- names(metacoefs$PM_Composition) <- 
  names(ses$PM_Composition) <- comp_names

# Compute RERs for IQR increase
rers <- Map(function(co, iqr) exp(co * iqr),
  metacoefs, iqrs)

# Compute confidence intervals
rercis <- Map(function(co, s, iqr) 
  exp(cbind(co - 1.96 * s, co + 1.96 * s) * iqr),
  metacoefs, ses, iqrs)

# Extract p-value
pvalues <- Map(function(co, s) 2 * (1 - pnorm(abs(co / s))),
  metacoefs, ses)

#----- Extract model criteria

# Extract Cochran's Q and I2
# Qs <- sapply(stage2res, function(x) summary(x)$qstat$Q)
# i2s <- sapply(stage2res, function(x) summary(x)$i2stat)

# Extract LRT p-value
lrts <- sapply(stage2res, function(x) lrt.mixmeta(x, stage2res$Null)$pvalue)
lrts["Null"] <- NA

# Extract information criterions
aics <- sapply(stage2res, AIC)
aiccs <- sapply(stage2res, AICc)
bics <- sapply(stage2res, BIC)

#----- Create final table

# Labels
modlabs <- c("Main", "Gas mixture", "Null", expression(O[x]), 
  expression(PM[2.5] ~ "Composition"))

# Bind all results into table
restab <- foreach(rer = rers, ci = rercis, p = pvalues, bic = bics, 
  lrt = lrts, aic = aiccs, lab = names(stage2res), .combine = bind_rows) %do% 
{
  # RER string
  rerstr <- if (nrow(ci) > 0) sprintf("%.4f (%.4f - %.4f)", 
    rer, ci[,1], ci[,2]) else NULL
  
  # Create table
  out <- cbind(model = lab, var = names(rer), rer = rerstr, 
    # q = q, i2 = i2, 
    lrt = lrt, aic = aic, bic = bic) |> as.data.frame()
}

# Column type
restab <- mutate(restab, across(all_of(c("lrt", "aic", "bic")), 
  as.numeric))

# Posterior probabilities from the BIC
deltabic <- exp(-(restab$bic - min(restab$bic)) / 2)
formatC(deltabic / sum(deltabic), format = "f")

#----- Output into a word file

# Labels
modlabs <- list(Ox = as_paragraph("O", as_sub("x")),
  PM_Composition = as_paragraph("PM", as_sub("2.5"), " Composition"),
  PMCI = as_paragraph("Main"), gas = as_paragraph("Gas Mixture"))
varlabs <- list(Ox = as_paragraph("O", as_sub("x")),
  SO4 = as_paragraph("SO", as_sub("4"), as_sup("2-")),
  NH4 = as_paragraph("NH", as_sub("4"), as_sup("+")),
  NIT = as_paragraph("NO", as_sub("3"), as_sup("-")),
  `log(I(PMCI + 1))` = as_paragraph("PMCI"),
  NO2_ppbv = as_paragraph("NO", as_sub("2")),
  SO2 = as_paragraph("SO", as_sub("2")),
  Ozone = as_paragraph("O", as_sub("3")),
  NH3 = as_paragraph("NH", as_sub("3"))
)

# Create flextable
resft <- flextable(restab[, 
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
for (m in seq_along(modlabs)) resft <- compose(resft, 
  j = "model", i = ~ model == names(modlabs)[m], value = modlabs[[m]])
for (v in seq_along(varlabs)) resft <- compose(resft, 
  j = "var", i = ~ var == names(varlabs)[v], value = varlabs[[v]])

# Add some horizontal borders
for (m in unique(restab$model)) resft <- hline(resft, 
  i = max(which(restab$model == m)),
  border = officer::fp_border(style = "dotted"))

resft <- resft |>
  # Relabel
  set_header_labels(model = "Model", var = "", #q = "Cochran's Q", 
    lrt = "LRT P-value", aic = "AIC", bic = "BIC", rer = "RER (95% CI)") |>
  # Border
  fix_border_issues() |>
  # Resize
  autofit() |> fit_to_width(20, unit = "cm")

save_as_docx(resft, path = "figures/Tab2_modelComparison.docx")

#-----------------------
# Association between RR and PMCI
#-----------------------

# Color palette
cntrpaldf <- subset(cities, conv) |>
  summarise(lon = mean(long), .by = country) |>
  left_join(unique(cities[, c("country", "countryname")])) |>
  mutate(isUSA = grepl("^USA", country)) |>
  group_by(isUSA) |>
  arrange(isUSA, lon) |>
  mutate(pal = scico(n(), palette = "batlow"))
cntrpal <- cntrpaldf$pal; names(cntrpal) <- cntrpaldf$countryname

# Shape palette
# shppal <- c(Canada = 16, UK = 17, Portugal = 10, Norway = 11,
#   Sweden = 15, Greece = 14, Romania = 18, China = 13)
shppal <- rep_len(15:18, sum(!cntrpaldf$isUSA))
names(shppal) <- subset(cntrpaldf, !isUSA, countryname, drop = T)
usashp <- rep_len(11:12, length(grep("^USA", cntrsum$countryname)))
names(usashp) <- subset(cntrpaldf, isUSA, countryname, drop = T)
shppal <- c(shppal, usashp)

# Plot association
rrplot <- ggplot(curvedata) + theme_bw() + 
  # geom_ribbon(aes(x = PMCI, ymin = ci.lb, ymax = ci.ub), alpha = .05) + 
  geom_line(aes(x = PMCI, y = fit), linewidth = 1) + 
  geom_line(aes(x = PMCI, y = ci.lb), linetype = 2) + 
  geom_line(aes(x = PMCI, y = ci.ub), linetype = 2) + 
  geom_point(aes(x = PMCI, y = blup, col = countryname, shape = countryname), 
    data = cities) + 
  geom_hline(yintercept = 1) + 
  scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
    name = "Country") + 
  scale_shape_manual(values = shppal, breaks = names(cntrpal), name = "Country") + 
  labs(x = "PMCI", y = "RR")

# Save
ggsave("figures/Fig1_RRpred.pdf", rrplot, height = 5, width = 8)



#--------------------------------
# Just for northern america
# 
# subset(cities, Region == "Northern America") |>
#   ggplot() + theme_bw() + 
#   geom_point(aes(x = PMCI, y = blup, col = countryname, shape = countryname)) + 
#   geom_hline(yintercept = 1) + 
#   scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
#     name = "Country") + 
#   scale_shape_manual(values = shppal, breaks = names(cntrpal), name = "Country") + 
#   labs(x = "PMCI", y = "RR (BLUP)")
# ggsave("temp/NA_RRBLUP.pdf")
# 
# 
# subset(cities, Region == "Northern America") |>
#   ggplot() + theme_bw() + 
#   geom_point(aes(x = PMCI, y = exp(coef), col = countryname, shape = countryname)) + 
#   geom_hline(yintercept = 1) + 
#   scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
#     name = "Country") + 
#   scale_shape_manual(values = shppal, breaks = names(cntrpal), name = "Country") + 
#   labs(x = "PMCI", y = "RR (1st stage)")
# ggsave("temp/NA_RR1st.pdf")
