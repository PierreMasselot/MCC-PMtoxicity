| Variable | Description |
|:---|:---|
|city | MCC city code |
| cityname | Name of the city |
| gdp | Country Gross Domestic Product |
| lat/long | Geographical coordinates |
| kgclzone/kgclzone1/kgclzone2 | levels of Köppen–Geiger climate classification |
| Region | regional classification of cities |
| country | MCC country code. Note that the US is further divided into regions |
| countryname | Name of the country. Note that the US is further divided into regions |
| PCA1-PCA4 | Principal components from the PCA on gaseous pollutants from Brook et al. PC1 is the basis for the PMCI. |
| NO2_molcm |  Average nitrogen dioxide concentration in molecules/cm2 |
| SO2 | Average sulfur dioxide concentration in Dobson unit |
| PM25 | Average fine particulate matter concentration in µg/m3 |
| Ozone | Average ozone concentration in ppbv |
| HCHO | Average formaldehyde concentration in molecules/cm2 |
| CO | Average carbon monoxide concentration in ppbv |
| NH3 | Average ammonia concentration in ppbv |
| NightLight | Night light |
| CAPI | Chronic air pollution index created in Brook et al. Corresponds to PC1 linearised between 0 and 100 which correspond to the pixels with lowest and highest PC1 values, respectively |
| PM25_Linear | PM25 linearised between 0 and 100 with the same method as the CAPI |
| PMCI | Pollutant Mixture Complexity Index |
| NO2_ppbv | Average nitrogen dioxide concentration in ppbv. Used to compute Ox. |
| Ox | Oxidative capacity of PM2.5 in ppbv |
| SO4 | Average proportion of PM2.5 that is sulfates |
| NH4 | Average proportion of PM2.5 that is ammonium |
| NIT | Average proportion of PM2.5 that is nitrates |
| BC | Average proportion of PM2.5 that is black carbon |
| OC | Average proportion of PM2.5 that is organic carbon |
| SS | Average proportion of PM2.5 that is sea salt |
| DUST | Average proportion of PM2.5 that is desert dust |
| GDPpc00 | City Gross Domestic Product in 2000 |
| GDPpc15 | City Gross Domestic Product in 2015 |
| E_GR_AV00 | NDVI in 2000 |
| E_GR_AV14 | NDVI in 2014 |
| B00 | Built-up environment in 2000 |
| B15 | Built-up environment in 2015 |
| tmean | Average temperature |
| trange | range of daily mean temperature |
| coef | PM2.5-related mortality log(RR) estimated in the first stage (quasi-Poisson regression) |
| v | variance of `coef` |
| conv | whether the quasi-Poisson regression fit converged |
| coefadj/vadj/convadj | same as above from first-stage models adjusted for O3 and NO2_molcm |
| deaths | city total deaths across the available period |
| periodmin | start of the studies period |
| periodmax | end of the studied period |
