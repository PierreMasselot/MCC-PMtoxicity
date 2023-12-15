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

#--------------------------
# Parameters
#--------------------------

#----- General

# Number of cores for parallelization
ncores <- max(1, detectCores() - 2)

#----- Data

# MCC countries
country_select <- c("aus8809", "can8615", "chi9615", "chl0414", "cyp0419", 
  "ecu1418", "est9718", "fnl9414", "fra0015", "ger9315", "grc0110", 
  "irn0215", "isr8520", "jap1115", "mex9814", "mlt9519", "nor6918", 
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

# ERF centering
cen <- 0

# Number of df per year
timedf <- 7

#----- Second stage analysis

# List of city-specific confounders
confvar <- c("GDPpc00", "GDPpc15", "E_GR_AV00", "E_GR_AV14", "B00", 
  "B15", "tmean", "trange", "PM2.5_Linear")
nconf <- length(confvar)

# Number of PCs
npc <- 2

# Ox weights
oxw <- c(no2 = 1.07, o3 = 2.075)
