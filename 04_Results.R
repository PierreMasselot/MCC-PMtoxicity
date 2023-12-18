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
  save_as_docx(path = "figures/countryDesc.docx")

#-----------------------
# Table with RER
#-----------------------

# Extract RER for each fitted model
nicetable <- foreach (fit = list(capifit, oxfit, adjfit), 
  lab = c("CAPI", "Ox", "CAPI (adjusted)"), 
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
  save_as_docx(path = "figures/TabRER.docx")

