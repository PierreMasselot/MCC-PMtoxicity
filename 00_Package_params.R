################################################################################
#
#                 MCC-PMToxicity
#             Part 0: Packages and parameters
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

#----- Analysis
library(tsModel)
library(dlnm)
library(mixmeta)
library(doParallel)
library(doSNOW)
library(splines)
library(compositions)
library(zCompositions)

#----- Output
library(knitr)
library(ggplot2)
library(flextable)
library(scico)
library(patchwork)
library(sf)

#--------------------------
# Parameters
#--------------------------

#----- General

# Number of cores for parallelization
ncores <- max(1, detectCores() - 2)

#----- Data

# MCC countries
country_select <- c("aus8819", "can8615", "chi9615", "chl0414", "cyp0419", 
  "ecu1418", "est9720", "fnl9414", "fra0017", "ger9320", "grc0110", 
  "irn0215", "isr8520", "jap1219", "mex9822", "mlt9519", "nor6918", "per0814",
  "por8018", "rom9416", "sa9713", "spa9014", "sui9513", "swe9010", "twn9414",
  "uk9020", "usa7306")

# Minimum record length available
nday <- 365

# Earliest year for analysis
yearstart <- 1999

#----- First-stage analysis

# Cross-basis parameters
maxlagp <- 1
arglagp <- list(fun = "strata") # Equivalent to MA
argvarp <- list(fun = "lin")

# Number of df per year
timedf <- 7

#----- Second stage analysis

# Formulas (mixed effect null model, and random effect)
nullform <- coef ~ pcs
ranform <- ~ 1|country/city

# Fit criterion
fitmethod <- "ml"

# List of city-specific confounders
confvar <- c("GDPpc00", "GDPpc15", "E_GR_AV00", "E_GR_AV14", "B00", 
  "B15", "tmean", "trange", "PM2.5_ug_m3")
nconf <- length(confvar)

# Number of PCs
npc <- 2

# Ox weights
oxw <- c(no2 = 1.07, o3 = 2.075)

# Function to perform Wald test
fwald <- function(full, null) {
  ind <- !names(coef(full)) %in% names(coef(null))
  coef <- coef(full)[ind]
  vcov <- vcov(full)[ind,ind]
  waldstat <- coef %*% solve(vcov) %*% coef
  df <- length(coef)
  pval <- 1 - pchisq(waldstat, df)
  return(list(waldstat = waldstat, pvalue = pval))
}

# Function to perform likelihood ratio test
lrt.mixmeta <- function(full, null) {
  llnull <- logLik(null)
  llfull <- logLik(full)
  lrstat <- -2 * (llnull - llfull)
  attributes(lrstat) <- NULL
  lrtdf <- attr(llfull, "df") - attr(llnull, "df")
  pval <- 1 - pchisq(lrstat, lrtdf)
  return(list(lrstat = lrstat, pvalue = pval))
}
