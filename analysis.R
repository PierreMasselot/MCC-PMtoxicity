################################################################################
#
#  MCC-PMToxicity project
#  
#  Air pollution mixture complexity and its effect on PM2.5-related mortality: 
#     a multi-country time-series study in 264 cities
#   
#  Environmental Epidemiology
#
#  MAIN SCRIPT
#
################################################################################

#--------------------------
# Packages
#--------------------------

#----- Data management
library(dplyr)
library(tibble)
library(stringr)
library(giscoR)
library(tidyr)

#----- Analysis
library(tsModel)
library(dlnm)
library(mixmeta)
library(doParallel)
library(doSNOW)
library(splines)
library(compositions)
library(zCompositions)
library(MuMIn) # For AICc

#----- Output
library(knitr)
library(ggplot2)
library(flextable)
library(scico)
library(patchwork)
library(sf)
library(ggnewscale)

#----- Source analysis functions
source("scripts/SecondStage_functions.R") # Functions to perform second-stage
source("scripts/Output_functions.R") # Tables and Figures
source("scripts/OtherRes_functions.R") # Additional results

#--------------------------
# Parameters
#--------------------------

#----- General

# Number of cores for parallelization
ncores <- max(1, detectCores() - 2)

#----- First-stage analysis

# Cross-basis parameters
maxlagp <- 1
arglagp <- list(fun = "strata") # Equivalent to MA
argvarp <- list(fun = "lin")

# Number of df per year
timedf <- 6

#----- Second stage analysis

# List of city-specific confounders
confvar <- c("GDPpc00", "GDPpc15", "E_GR_AV00", "E_GR_AV14", "B00", 
  "B15", "tmean", "trange", "PM25")

# Number of PCs and names
npc <- 2
pcnms <- sprintf("conf_pc%i", 1:npc)

# Formulas (mixed effect null model, and random effect)
nullform <- sprintf("coef ~ %s", paste(pcnms, collapse = " + ")) |>
  as.formula()
ranform <- ~ 1|country/city

# Fitting of the second-stage meta-regression
fitmethod <- "ml"

#----- Outputing

# Model labels for the output table
modlabs <- list(Ox = as_paragraph("O", as_sub("x")),
  PMcomposition = as_paragraph("PM", as_sub("2.5"), " Composition"),
  PMCI = as_paragraph("Main"), Gas = as_paragraph("Gas Mixture"))

# Variable labels for the output table
varlabs <- list(
  Ox = as_paragraph("O", as_sub("x")),
  SO4 = as_paragraph("SO", as_sub("4"), as_sup("2-")),
  NH4 = as_paragraph("NH", as_sub("4"), as_sup("+")),
  NIT = as_paragraph("NO", as_sub("3"), as_sup("-")),
  `log(I(PMCI + 1))` = as_paragraph("PMCI"),
  NO2_ppbv = as_paragraph("NO", as_sub("2")),
  SO2 = as_paragraph("SO", as_sub("2")),
  Ozone = as_paragraph("O", as_sub("3")),
  NH3 = as_paragraph("NH", as_sub("3"))
)


###################################################
# Main analysis
###################################################

#--------------------------
# First-stage analysis
#--------------------------

# Due to restricted data-sharing agreements, the first-stage is non-reproducible
# The result of the first-stage is saved in folder "data"

# This section checks if original data are there and if not, 
#   load the pre-saved first-stage results
if(dir.exists("original")){
  
  # Link MCC and the various pollution data together
  source("original/LinkData.R")
  
  # Run the first-stage analysis (the code is available)
  source("scripts/FirstStage_script.R")

} else {
  # If not available, load data
  cities <- read.csv("data/citydata.csv")
}

#--------------------------
# Second-stage meta-regression and results
# This part is fully reproducible
#--------------------------

#----- Main model

# Compute the PCs on confounders and add to data
pcs <- run_pca(cities, 
  confvar = c("GDPpc00", "GDPpc15", "E_GR_AV00", "E_GR_AV14", "B00", 
    "B15", "tmean", "trange", "PM25"), 
  npc = npc)
colnames(pcs) <- pcnms
cities <- cbind(cities, pcs)


# Run the second-stage meta-regression model
mainmod <- run_metareg(cities, nullform = nullform, 
  mainvar = ~ log(I(PMCI + 1)), random = ranform, S = cities$v, 
  method = fitmethod, subset = cities$conv)

#----- Model comparison

# List of models defined by formula
formlist <- list(
  Null = ~ 1, 
  Gas = ~ NO2_ppbv + SO2 + Ozone + HCHO + CO + NH3,
  Ox = ~ Ox)

# Run alternative models
othermod <- lapply(formlist, run_metareg, cities = cities, nullform = nullform,
  random = ranform, S = cities$v, method = fitmethod, subset = cities$conv)

# Run composition model (slightly different extraction of results)
othermod$PMcomposition <- run_metareg_comp(
  mainvar = ~ alr(cbind(SO4, NH4, NIT, BC, OC, SS, DUST)), 
  cities = cities, nullform = nullform, random = ranform,
  S = cities$v, method = fitmethod, subset = cities$conv)

#--------------------------
# Tables and Plots
#--------------------------

#----- Table 1: country descriptive statistics
table_descriptive(subset(cities, conv), "figures/Tab1_countryDesc.docx")

#----- Table 2: model comparison
table_comparison(c(list(Main = mainmod), othermod), 
  modlabs = modlabs, varlabs = varlabs, 
  cols = c("model", "var", "rer", "LRT", "AICc", "BIC"),
  file = "figures/Tab2_modelComparison.docx")

#----- Figure 1: RR vs PMCI

# Create color and point shape palettes
pals <- make_palettes(cities)

# Plot
plot_RRpred(list(Main = mainmod), var = "PMCI", cities = cities, 
  countrypal = pals, modelpal = "black", 
  file = "figures/Fig1_RRpred.pdf")

###################################################
# Sensitivity analysis: NO2 and O3 adjusted RR
###################################################

#--------------------------
# Rerun meta-regression models with adjusted RRs and variances
#--------------------------

# Formulas
sensform <- update(nullform, coefadj ~ .)
formlist2 <- c(list(Main = ~ log(I(PMCI + 1))), formlist)

# Fit models (PM composition using a different function to extract results)
sensmod <- lapply(formlist2, run_metareg, cities = cities, nullform = sensform,
  random = ranform, S = cities$vadj, method = fitmethod, 
  subset = cities$convadj)
sensmod$PMcomposition <- run_metareg_comp(
  mainvar = ~ alr(cbind(SO4, NH4, NIT, BC, OC, SS, DUST)), 
  cities = cities, nullform = sensform, random = ranform,
  S = cities$vadj, method = fitmethod, subset = cities$convadj)

#--------------------------
# Tables and Plots
#--------------------------

# Description of data
table_gasDesc(subset(cities, convadj), 
  "figures/SupTab1_NO2O3_countryDesc.docx")

# Comparison of RRs
figure_rrComparison(cities, palettes = pals, file = "figures/SupFig4_RRadj.png")

# Comparison of models
table_comparison(sensmod, modlabs = modlabs, varlabs = varlabs, 
  cols = c("model", "var", "rer", "LRT", "AICc", "BIC"),
  file = "figures/SupTab5_sensitivityAnalysis.docx")

###################################################
# Additional results
###################################################

#----- Maps of pollutant levels

figure_polmap(
  vars = c("PMCI", "PM25", "NO2_ppbv", "Ozone", "SO2", "HCHO", "CO", "NH3"),
  varlabs = c("PMCI", expression(PM[2.5] ~ "(" * mu * "g/m"^{3} * ")"),
    expression(NO[2] ~ "(ppbv)"), "Ozone (ppbv)", SO[2] ~ "(du)", 
    expression("HCHO (mol / cm"^{2} * ")"), "CO (ppbv)", 
    expression(NH[3] ~ "(ppbv)")),
  transvars = c("identity", "log10", "log10", "identity", "log10", "log10", 
    "log10", "log10"),
  file = "figures/SupFig1_pollmaps.png", height = 20
)

#----- Selection of main model

# List of formulas
pmciterms <- list(
  Linear = ~ PMCI, 
  Log = ~ log(I(PMCI + 1)), 
  `ns (df = 2)` = ~ ns(PMCI, df = 2),
  `ns (df = 3)` = ~ ns(PMCI, df = 3)
)

# Fit models
selres <- lapply(pmciterms, run_metareg, cities = cities, nullform = nullform,
  random = ranform, S = cities$v, method = fitmethod, subset = cities$conv)

# Comparison table
table_comparison(selres, modlabs = NULL, varlabs = NULL, 
  cols = c("model", "AICc", "BIC", "psi_country", "psi_city"),
  file = "figures/SupTab3_modelComparison.docx")

# Plot predicted RR
modpal <- scico(length(selres), palette = "glasgow")
names(modpal) <- names(selres)
plot_RRpred(selres, var = "PMCI", cities = cities, modelpal = modpal,
  file = "figures/SupFig2_ComparisonCurve.png")

#----- Second-stage meta-regression heterogeneity

# Create table as in main manuscript
table_comparison(c(list(Main = mainmod), othermod), 
  modlabs = modlabs, cols = c("model", "Q", "i2", "psi_country", "psi_city"),
  file = "figures/SupTab4_metaCrit.docx")

#----- Residuals

# Extract residuals
cities[cities$conv, "resid"] <- residuals(mainmod$fit)

# Plot
figure_residuals(cities, residvar = "resid",
  xvars = c("PMCI", "countryname", "deaths", "PM25"),
  xlabs = c("PMCI", "", "Deaths", expression(PM[2.5]~"("*mu*g*"/"*m^{3}*")")),
  smooth = c(T, F, T, T), countrypal = pals, log = c(F, F, T, T),
  file = "figures/SupFig3_Residuals.png", width = 10, height = 10)

