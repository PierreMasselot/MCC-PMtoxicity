#####################################################################
#
#                 MCC HetPoll
#               Part 4: Plots
#
#####################################################################

#-----------------------
# Description of cities
#-----------------------

# Number of countries
length(unique(cities$country))

# Extract number of cities by country
table(subset(cities, conv & !is.na(ox), countryname, drop = T)) |> 
  as.data.frame() |>
  flextable() |>
  delete_part() |>
  autofit() |>
  save_as_docx(path = "figures/countryDesc_sameDat.docx")

#-----------------------
# Table with RER
#-----------------------

# Extract RER for each fitted model
nicetable <- foreach (fit = list(capifit, oxfit), lab = c("CAPI", "Ox"), 
  .combine = rbind) %do% 
{
  
  # Extract results and convert to RER
  coefres <- summary(fit)$coefficients[2, 
    c("Estimate", "95%ci.lb", "95%ci.ub", "Pr(>|z|)"), drop = F] |> 
    as.data.frame() |>
    rename("CI_lower" = "95%ci.lb", "CI_upper" = "95%ci.ub", 
      "P-value" = "Pr(>|z|)") |>
    mutate(
      Estimate = exp(Estimate * IQR(model.matrix(fit)[,2])), 
      CI_lower = exp(CI_lower * IQR(model.matrix(fit)[,2])),
      CI_upper = exp(CI_upper * IQR(model.matrix(fit)[,2])))
  
  # Format numbers
  coefres <- mutate(coefres, variable = lab,
    RER = sprintf("%1.4f (%1.4f to %1.4f)", Estimate, CI_lower, CI_upper),
    `P-value` = sprintf("%1.4f", `P-value`)) 
}

# Output into a word file
flextable(nicetable[,c("variable", "RER", "P-value")]) |>
  set_header_labels(variable = "") |>
  autofit() |>
  save_as_docx(path = "figures/TabRER_sameDat.docx")

#----------------------------
# Plot how RR affects CAPI
#----------------------------

#----- Predictions for various values

# Prediction data.frame
capi_pred <- cities[, c(attr(terms(capi), "term.labels"), "coef")] |>
  mutate(PM2.5_Linear = mean(PM2.5_Linear), gdp = mean(gdp))

# Prediction
capi_RR <- predict(capi, capi_pred, se = T) |>
  cbind(capi_pred)

# CIs and RR
capi_RR <- mutate(capi_RR, RR = exp(fit), low = exp(fit - 1.96 * se), 
  high = exp(fit + 1.96 * se))

#----- Plot

ggplot(capi_RR, aes(x = PM2.5_toxicity)) + theme_bw() + 
  geom_point(aes(y = exp(coef))) + 
  geom_line(aes(y = RR)) + 
  geom_line(aes(y = low), linetype = 2) + 
  geom_line(aes(y = high), linetype = 2) + 
  geom_hline(yintercept = 1) + 
  labs(y = "RR")

ggsave("figures/RR_v_toxi.pdf")
