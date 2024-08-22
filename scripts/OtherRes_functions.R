#####################################################################
#
#                 MCC-PMToxicity
#            Other results
#
#####################################################################

##############################################
# Used in the sensitivity analysis
##############################################

#--------------------------
# NO2 and O3 description
#--------------------------

table_gasDesc <- function(cities, file = NULL) {

  # Compute summaries of pollutants
  senscntrsum <- cities |>
    summarise(ncities = n(), no2m = mean(NO2_ppbv), 
      no2q1 = quantile(NO2_ppbv, .25), no2q3 = quantile(NO2_ppbv, .75),
      o3m = mean(Ozone), o3q1 = quantile(Ozone, .25), 
      o3q3 = quantile(Ozone, .75),
      .by = countryname) |>
    arrange(countryname)
  
  # Add total
  senstotsum <- cities |>
    summarise(ncities = n(), no2m = mean(NO2_ppbv), 
      no2q1 = quantile(NO2_ppbv, .25), no2q3 = quantile(NO2_ppbv, .75),
      o3m = mean(Ozone), o3q1 = quantile(Ozone, .25), 
      o3q3 = quantile(Ozone, .75))
  senscntrsum <- bind_rows(senscntrsum, senstotsum)
  
  # Add name
  senscntrsum <- merge(senscntrsum, 
    subset(countries, select = c(country, countryname)),
    all.x = T)
  senscntrsum$countryname[is.na(senscntrsum$countryname)] <- "Total"
  
  # Create pretty column
  senscntrsum <- mutate(senscntrsum, 
    NO2 = sprintf("%.2f (%.2f - %.2f)", no2m, no2q1, no2q3),
    O3 = sprintf("%.2f (%.2f - %.2f)", o3m, o3q1, o3q3)) |>
    subset(select = c(countryname, ncities, NO2, O3))
  
  #----- Create and export word table
  senssumtab <- flextable(senscntrsum) |>
    set_header_labels(countryname = "Country", ncities = "Number of cities") |>
    compose(part = "header", j = "NO2", 
      value = as_paragraph("Average NO", as_sub("2"), " in ppbv (IQR)")) |>
    compose(part = "header", j = "O3", 
      value = as_paragraph("Average O", as_sub("3"), " in ppbv (IQR)")) |>
    bold(i = nrow(senscntrsum)) |>
    autofit() |> width(j = 2, width = .75)
  
  # Save and return
  if(!is.null(file)) save_as_docx(senssumtab, path = file)

  return(invisible(senssumtab))
}


#--------------------------
# Comparison of RRs when adjusted by NO2/O3
#--------------------------

figure_rrComparison <- function(cities, palettes, file = NULL) {
  
  # List of featured countries
  redcntr <- subset(cities, !is.na(coefadj), countryname, drop = T) |> unique()
  senspal <- palettes$color[names(palettes$color) %in% redcntr]
  sensshp <- palettes$shape[names(palettes$shape) %in% redcntr]
  
  # Number of lower RRs
  summarise(cities, lower = mean(coefadj < coef, na.rm = T), 
    meandiff = mean(exp(coef) - exp(coefadj), na.rm = T))
  mutate(cities, diff = abs(exp(coef) - exp(coefadj))) |>
    arrange(desc(diff)) |>
    head()
  
  # Plot
  rrcompplot <- ggplot(cities) + theme_bw() + 
    geom_point(aes(x = exp(coef), y = exp(coefadj), 
      col = countryname, shape = countryname)) + 
    geom_hline(yintercept = 1) + 
    geom_vline(xintercept = 1) + 
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = senspal, breaks = names(senspal), 
      name = "Country", drop = T) + 
    scale_shape_manual(values = sensshp, breaks = names(senspal), 
      name = "Country", drop = T) + 
    labs(x = "RR (main)", y = "RR (adjusted)")
  
  if (!is.null(file)) ggsave(file, rrcompplot, height = 5, width = 8)
  
  return(invisible(rrcompplot))
}

##############################################
# Other results
##############################################

#--------------------------
# Maps of pollutant levels
#--------------------------

figure_polmap <- function(vars, varlabs, transvars, file = NULL, ...){
  
  # Load world boundaries
  cntrbnd <- gisco_get_countries(resolution = "20")
  
  # Loop
  mapplots <- foreach(v = vars, lab = varlabs, tr = transvars) %do% {
    ggplot(cities) + theme_void() +
      theme(panel.grid.major = element_line(color = "lightgrey")) +
      geom_sf(data = cntrbnd, fill = grey(.95), col = grey(.5)) + 
      geom_point(aes(x = long, y = lat, col = .data[[v]]), alpha = .5, size = 2) +
      coord_sf(crs = "+proj=robin", default_crs = 4326, 
        ylim = c(0, 90)) + 
      scale_color_viridis_c(name = "", transform = tr) +
      labs(title = lab)
  }
  
  # Put together
  pollmaps <- wrap_plots(mapplots, ncol = 1)
  
  # Save and return
  if(!is.null(file)) ggsave(file, pollmaps, ...)
  
  return(invisible(pollmaps))
}

#--------------------------
# Plotting residuals
#--------------------------

figure_residuals <- function(cities, residvar, xvars, xlabs, smooth = T, 
  countrypal, log = F, file = NULL, ...)
{
  
  # Create plots for all X variables
  allplots <- foreach(var = xvars, lab = xlabs, 
    sm = rep_len(smooth, length(xvars)), l = log) %do% 
  {
    g <- ggplot(cities) + theme_bw()
    if (sm) g <- g + geom_smooth(aes(x = .data[[var]], y = .data[[residvar]]), 
      method = "gam", col = 1)
    g <- g + 
      geom_point(aes(x = .data[[var]], y = .data[[residvar]], col = countryname, 
        shape = countryname)) + 
      geom_hline(yintercept = 0) +
      scale_color_manual(values = countrypal$color, 
        breaks = names(countrypal$color), name = "Country") + 
      scale_shape_manual(values = countrypal$shape, 
        breaks = names(countrypal$shape), name = "Country") + 
      labs(x = lab, y = "log(RR) residual")
    if (l) g <- g + scale_x_log10()
    if (!is.numeric(cities[[var]])) g <- g + 
      theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))
    g
  }
  
  # Put plots together
  gall <- wrap_plots(allplots, ncol = 2, guides = "collect") & 
    theme(legend.position = "bottom", legend.title.position = "top")
  
  # Save and return 
  if(!is.null(file)) ggsave(file, gall, ...)
  
  return(invisible(gall))
}
