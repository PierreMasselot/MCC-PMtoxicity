#####################################################################
#
#                 MCC HetPoll
#               Part 7: Supplementary plots
#
#####################################################################

#----------------------------
# Pollutant mix description
#----------------------------

#----- World maps

# Variables to plot
pmcivars <- c("PMCI", "PM25", "NO2_ppbv", "Ozone", "SO2", "HCHO", "CO", "NH3")
varlabs <- c("PMCI", expression(PM[2.5] ~ "(" * mu * "g/m"^{3} * ")"),
  expression(NO[2] ~ "(ppbv)"), "Ozone (ppbv)", SO[2] ~ "(du)", 
  expression("HCHO (mol / cm"^{2} * ")"), "CO (ppbv)", 
  expression(NH[3] ~ "(ppbv)"))
transvars <- c("identity", "log10", "log10", "identity", "log10", "log10", 
  "log10", "log10")

# Load world boundaries
cntrbnd <- gisco_get_countries(resolution = "20")

# Loop
mapplots <- foreach(v = pmcivars, lab = varlabs, tr = transvars) %do% {
  ggplot(cities) + theme_void() +
    theme(panel.grid.major = element_line(color = "lightgrey")) +
    geom_sf(data = cntrbnd, fill = grey(.95), col = grey(.5)) + 
    geom_point(aes(x = long, y = lat, col = .data[[v]]), alpha = .5, size = 2) +
    coord_sf(crs = "+proj=robin", default_crs = 4326, 
      ylim = c(0, 90)) + 
    scale_color_viridis_c(name = "", transform = tr) +
    labs(title = lab)
}

# Put together and save
pollmaps <- wrap_plots(mapplots, ncol = 1)
ggsave("figures/SupFig1_pollmaps.png", pollmaps, height = 20)

#----------------------------
# Model selection
#----------------------------

#----- Compare the models

# Terms for PMCI
modterms <- list(Linear = ~ . + PMCI, 
  Log = ~ . + log(I(PMCI + 1)), 
  `ns (df = 2)` = ~ . + ns(PMCI, df = 2),
  # `ns (df = 2)` = ~ . + ns(PMCI, knots = 1),
  `ns (df = 3)` = ~ . + ns(PMCI, df = 3)
  # `ns (k = 1,2)` = ~ . + ns(PMCI, knots = c(1, 2))
)

# Formulas
modforms <- lapply(modterms, update, old = nullform)

# Apply models
modtests <- lapply(modforms, function(f) mixmeta(f, S = v, random = ranform, 
  data = cities, method = fitmethod, subset = conv))

# Extract result (random effect and information criteria)
modres <- foreach(mod = modtests, lab = names(modterms), 
  .combine = rbind) %do% 
{
  cbind(data.frame(Model = lab, AIC = AICc(mod), BIC = BIC(mod)), 
    sqrt(t(unlist(mod$Psi))))
}

# Export as a table
comparft <- flextable(modres) |>
  set_header_labels(country = "Country", city = "City") |>
  add_header_row(values = c("", "Random effect std. deviation"),
    colwidths = c(3, 2), top = TRUE) |>
  align(align = "center", part = "header", j = 1:5) |>
  colformat_double(j = c("BIC", "AIC"), digits = 2, big.mark = "") |>
  colformat_double(j = c(4, 5), digits = 6) |>
  autofit()
save_as_docx(comparft, path = "figures/SupTab3_modelComparison.docx")
  
# Display posterior probabilities from BIC
deltabic <- exp(-(modres$BIC - min(modres$BIC)) / 2)
deltabic / sum(deltabic)


#----- Check the predicted curves

# Make predictions
predtests <- lapply(modtests, predict, curvedata, ci = T) |> 
  do.call(what = cbind)

# Color palette
predpal <-  scico(length(modtests), palette = "glasgow")
  
# Plot predictions
matplot(curvedata$PMCI, exp(predtests), type = "l", 
  col = rep(predpal, each = 3),
  lty = rep(c(1, 3), 1:2), lwd = rep(c(2, 1), 1:2),
  xlab = "PMCI", ylab = "RR")
abline(h = 1)
legend("bottomright", legend = names(modtests), col = predpal, lty = 1, lwd = 2,
  bty = "n", ncol = 2)

dev.print(png, filename = "figures/SupFig2_ComparisonCurve.png", units = "in",
  res = 300)

#----- Check residuals

# # Extract residuals
# residtests <- lapply(modtests, residuals, curvedata, ci = T)
# 
# # Plot residuals
# par(mfrow = n2mfrow(length(residtests), asp = 2/3))
# Map(function(r, lab) {
#   plot(r ~ PMCI, data = subset(cities, conv), pch = 16, 
#     ylab = "log(RR) residuals", main = lab)
#   abline(h = 0)
# }, residtests, names(residtests))
# 
# dev.print(png, filename = "figures/SupFig3_ComparisonResid.png", units = "in",
#   res = 300)

#----------------------------
# Additional results
#----------------------------

# Extract supplementary results
supres <- foreach(mod = stage2res, .combine = rbind) %do% {
  data.frame(Q = summary(mod)$qstat$Q, i2 = summary(mod)$i2stat,
    lapply(mod$Psi, sqrt))
}
names(supres)[3:4] <- c("country", "city")
supres$mod <- names(stage2res)

# Export
supft <- flextable(supres[, c("mod", "Q", "i2", "country", "city")]) |>
  colformat_double(j = c("Q", "i2"), digits = 2, big.mark = "") |>
  colformat_double(j = c("country", "city"), digits = 4, big.mark = "") |>
  align(j = "mod", align = "right") |>
  compose(j = "i2", value = as_paragraph("I", as_sup("2")), part = "header")
for (m in seq_along(modlabs)) supft <- compose(supft, 
  j = "mod", i = ~ mod == names(modlabs)[m], value = modlabs[[m]])
supft <- set_header_labels(supft, mod = "Model", Q = "Cochran's Q", 
    country = "Country Std. Dev.", city = "City Std. Dev.") |>
  autofit()
save_as_docx(supft, path = "figures/SupTab4_metaCrit.docx")

#----------------------------
# Residuals
#----------------------------

## Reuses palettes from Figure 1

# Plot residuals vs PMCI
residpmci <- ggplot(cities) + theme_bw() + 
  geom_smooth(aes(x = PMCI, y = mar_resid), method = "gam", col = 1) + 
  geom_point(aes(x = PMCI, y = mar_resid, col = countryname, 
    shape = countryname)) + 
  geom_hline(yintercept = 0) +
  scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
    name = "Country") + 
  scale_shape_manual(values = shppal, breaks = names(shppal), 
    name = "Country") + 
  labs(x = "PMCI", y = "log(RR) residual")

# Plot residuals vs country
residcntr <- ggplot(cities) + theme_bw() + 
  geom_point(aes(x = factor(countryname, levels = cntrpaldf$countryname), 
    y = mar_resid, col = countryname, shape = countryname)) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
    name = "Country") + 
  scale_shape_manual(values = shppal, breaks = names(shppal), 
    name = "Country") + 
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(x = "", y = "log(RR) residual")

# Plot residuals vs total mortality
residdeaths <- ggplot(cities) + theme_bw() + 
  geom_smooth(aes(x = deaths, y = mar_resid), method = "gam", col = 1) + 
  geom_point(aes(x = deaths, y = mar_resid, col = countryname, 
    shape = countryname)) + 
  geom_hline(yintercept = 0) +
  scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
    name = "Country") + 
  scale_shape_manual(values = shppal, breaks = names(shppal), 
    name = "Country") + 
  scale_x_log10() +
  labs(x = "Deaths", y = "log(RR) residual")

# Plot residuals vs PM2.5 level
residpm25 <- ggplot(cities) + theme_bw() + 
  geom_smooth(aes(x = PM25, y = mar_resid), method = "gam", col = 1) + 
  geom_point(aes(x = PM25, y = mar_resid, col = countryname, 
    shape = countryname)) + 
  geom_hline(yintercept = 0) +
  scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
    name = "Country") + 
  scale_shape_manual(values = shppal, breaks = names(shppal), 
    name = "Country") + 
  scale_x_log10() +
  labs(x = expression(PM[2.5]~"("*mu*g*"/"*m^{3}*")"), y = "log(RR) residual")

# Put together
(residpmci + residpm25) / (residdeaths + residcntr) + 
  plot_layout(guides = "collect")

# Save
ggsave("figures/SupFig3_Residuals.png")
