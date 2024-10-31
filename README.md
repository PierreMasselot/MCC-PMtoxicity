# MCC-PMtoxicity

R code and intermediary data performing the analysis from the paper **Air pollution mixture complexity and its effect on PM2.5-related mortality: a multi-country time-series study in 264 cities** published in *Environmental Epidemiology* (https://doi.org/10.1097/EE9.0000000000000342). 

## Data

Due to a restrictive data sharing agreement between MCC members on daily mortality and air pollutants, the original data cannot be shared publicly. We instead share city-level aggregated data as well as the log relative risks of PM2.5-related mortality estimated in the first-stage of the analysis. The city-level data can be found in the `data` folder, along with a codebook describing all variables.

## Analysis

The file `analysis.R` is the master file running the full analysis. It calls functions from the R files in the `scripts` folder. 

Due to the restrictions mentioned above, the first-stage of the analysis is not reproducible. Instead, we saved the estimated city-specific log(RR) along with city-specific characteristic and any useful information in the `data` folder. The second-stage of the analysis is thus reproducible using the saved first-stage results. The script `analysis.R` automatically skips the first-stage of the analysis and load the saved results in `data`, unless the use has access to the original data, which would be stored in a `original` folder.

Note that the script used to perform the first stage is available and can therefore be checked and copied by any interested reader of the study.