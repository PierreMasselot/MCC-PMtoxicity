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
  summarise(ncities = n(), PM25 = mean(PM2.5_ug_m3), 
    PM25q1 = quantile(PM2.5_ug_m3, .25), PM25q3 = quantile(PM2.5_ug_m3, .75),
    pmci = mean(PMCI), pmciq1 = quantile(PMCI, .25), 
    pmciq3 = quantile(PMCI, .75),
    deaths = sum(deaths), 
    periodmin = min(periodmin), periodmax = max(periodmax),
  .by = countryname) |>
  arrange(countryname)

# Add total
totsum <- subset(cities, conv) |>
  summarise(ncities = n(), PM25 = mean(PM2.5_ug_m3), 
    PM25q1 = quantile(PM2.5_ug_m3, .25), PM25q3 = quantile(PM2.5_ug_m3, .75),
    pmci = mean(PMCI),
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
    PM25 = sprintf("%.2f (%.2f - %.2f)", PM25, PM25q1, PM25q3),
    PMCI = sprintf("%.2f (%.2f - %.2f)", pmci, pmciq1, pmciq3),
    period = sprintf("%i - %i", periodmin, periodmax)) |>
  subset(select = c(countryname, ncities, period, deaths, PM25, PMCI))

#----- Create and export word table
sumtab <- flextable(cntrsum) |>
  set_header_labels(countryname = "Country", ncities = "Number of cities",
    period = "Period", deaths = "Total mortality", 
    PMCI = "Average PMCI (IQR)") |>
  compose(part = "header", j = 5, 
    value = as_paragraph("Average PM", as_sub("2.5"), " (IQR)")) |>
  bold(i = nrow(cntrsum)) |>
  autofit()
save_as_docx(sumtab, path = "figures/Tab1_countryDesc.docx")

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
Qs <- sapply(stage2res, function(x) summary(x)$qstat$Q)
i2s <- sapply(stage2res, function(x) summary(x)$i2stat)

# Extract LRT p-value
lrts <- sapply(stage2res, function(x) lrt.mixmeta(x, stage2res$Null)$pvalue)

# Extract AICs
aics <- sapply(stage2res, AIC)

#----- Create final table

# Labels
modlabs <- c("Main", "Null", "Ox", "PM Composition")

# Bind all results into table
restab <- foreach(rer = rers, ci = rercis, p = pvalues, q = Qs, i2 = i2s, 
  lrt = lrts, aic = aics, lab = modlabs, .combine = bind_rows) %do% 
{
  # RER string
  rerstr <- if (nrow(ci) > 0) sprintf("%.4f (%.4f - %.4f)", 
    rer, ci[,1], ci[,2]) else NULL
  
  # Create table
  out <- cbind(model = lab, var = names(rer), rer = rerstr, 
    q = q, i2 = i2, lrt = lrt, aic = aic) |> as.data.frame()
  
}

# Column type
restab <- mutate(restab, across(all_of(c("q", "i2", "lrt", "aic")), 
  as.numeric))

# Change some text
restab[restab$model == "Main", "var"] <- "PMCI"

#----- Output into a word file
resft <- flextable(restab[, 
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
  colformat_double(j = c("q", "i2", "aic"), digits = 2, big.mark = "") |>
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
ggsave("figures/Fig1_RRpred.pdf", rrplot, height = 5)
