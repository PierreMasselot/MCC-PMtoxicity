#####################################################################
#
#                 MCC-PMToxicity
#            Plotting functions
#
#####################################################################

#--------------------------
# Descriptive table
#--------------------------

table_descriptive <- function(cities, file = NULL) {
  
  #----- Summary by country
  
  # Compute summaries of pollutants
  cntrsum <- cities |>
    summarise(ncities = n(), PM25m = mean(PM25), 
      PM25q1 = quantile(PM25, .25), PM25q3 = quantile(PM25, .75),
      pmcim = mean(PMCI), pmciq1 = quantile(PMCI, .25), 
      pmciq3 = quantile(PMCI, .75),
      deaths = sum(deaths), 
      periodmin = min(periodmin), periodmax = max(periodmax),
      .by = countryname) |>
    arrange(countryname)
  
  # Add total
  totsum <- cities |>
    summarise(ncities = n(), PM25m = mean(PM25), 
      PM25q1 = quantile(PM25, .25), PM25q3 = quantile(PM25, .75),
      pmcim = mean(PMCI),
      pmciq1 = quantile(PMCI, .25), 
      pmciq3 = quantile(PMCI, .75),
      deaths = sum(deaths), 
      periodmin = min(periodmin), periodmax = max(periodmax))
  cntrsum <- bind_rows(cntrsum, totsum)
  
  # Add name
  cntrsum$countryname[is.na(cntrsum$countryname)] <- "Total"
  
  # Create pretty column
  cntrsum <- mutate(cntrsum, 
    PM25 = sprintf("%.2f (%.2f - %.2f)", PM25m, PM25q1, PM25q3),
    PMCI = sprintf("%.2f (%.2f - %.2f)", pmcim, pmciq1, pmciq3),
    period = sprintf("%i - %i", periodmin, periodmax)) |>
    subset(select = c(countryname, ncities, period, deaths, PM25, PMCI))
  
  #----- Create and export word table
  sumtab <- flextable(cntrsum) |>
    set_header_labels(countryname = "Country", ncities = "Number of cities",
      period = "Period", deaths = "Total mortality", 
      PMCI = "Average PMCI (IQR)") |>
    compose(part = "header", j = 5, 
      value = as_paragraph("Average PM", as_sub("2.5"), " in µg/m", 
        as_sup("3"), " (IQR)")) |>
    bold(i = nrow(cntrsum)) |>
    autofit() |> 
    fit_to_width(20, unit = "cm")
  
  # Save and return table
  if (!is.null(file)) save_as_docx(sumtab, path = file)
  
  return(invisible(sumtab))
}


#--------------------------
# Model comparison table
#--------------------------

table_comparison <- function(modlist, modlabs, varlabs, cols = colnames(restab), 
  file = NULL)
{
  
  if (!"var" %in% cols) varlabs <- NULL
  
  #----- Put results into a table
  
  restab <- foreach(mod = modlist, lab = names(modlist),
    .combine = bind_rows) %do% 
  {
    
    # RER with CI
    rerstr <- if (nrow(mod$RER) > 0) sprintf("%.4f (%.4f - %.4f)", 
      mod$RER[, "rer"], mod$RER[,"low"], mod$RER[,"high"]) else NULL
    
    # Random effect variances
    psi <- sqrt(t(unlist(mod$fit$Psi)))
    colnames(psi) <- sprintf("psi_%s", colnames(psi))
    
    # Put into a table
    bind_cols(model = lab, var = rownames(mod$RER), rer = rerstr, 
        t(mod$score[c("LRT", "AICc", "BIC")]), 
        Q = summary(mod$fit)$qstat$Q, i2 = summary(mod$fit)$i2stat,
        psi) |> 
      as.data.frame()
  }
  
  # # Column type
  # restab <- mutate(restab, across(all_of(c("LRT", "AICc", "BIC")), 
  #   as.numeric))
  
  # Select columns and filter rows if necessary
  restab <- restab[, cols]
  restab <- restab[!duplicated(restab),]
  
  #----- Output into a word file
  
  # Create flextable
  resft <- flextable(restab[, cols]) |>
    
    # Merge identical rows
    merge_v() |>
    valign(valign = "top")
    
  # Put the best value in bold
  for (j in intersect(cols, c("LRT", "AICc", "BIC"))) resft <- bold(resft,
    as.formula(sprintf("~ %s == min(%s, na.rm = T)", j, j)), j)
    
  # Format columns
  resft <- resft |>
    colformat_double(j = intersect(cols, c("AICc", "BIC", "Q", "i2")), 
      digits = 2, big.mark = "") |>
    colformat_double(j = intersect(cols, c("LRT", "psi_country", "psi_city")), 
      digits = 4, na_str = "-")
  
  # Change some labels
  for (m in seq_along(modlabs)) resft <- compose(resft, 
    j = "model", i = as.formula(sprintf("~ model == '%s'", names(modlabs)[m])), 
    value = modlabs[[m]])
  for (v in seq_along(varlabs)) resft <- compose(resft, 
    j = "var", i = as.formula(sprintf("~ var == '%s'", names(varlabs)[v])), 
    value = varlabs[[v]])
  
  # Add some horizontal borders
  for (m in unique(restab$model)) resft <- hline(resft, 
    i = max(which(restab$model == m)),
    border = officer::fp_border(style = "dotted"))
  
  # Relabel
  if ("i2" %in% cols) resft <- compose(resft, j = "i2", 
    value = as_paragraph("I", as_sup("2")), part = "header") 
  resft <- set_header_labels(resft, model = "Model", var = "", 
    LRT = "LRT P-value", AICc = "AIC", BIC = "BIC", rer = "RER (95% CI)",
    Q = "Cochran's Q", psi_country = "Country Std. Dev.", 
    psi_city = "City Std. Dev.")
    
  # Final layout
  resft <- resft |>
    # Border
    fix_border_issues() |>
    
    # Resize
    autofit() |> fit_to_width(20, unit = "cm")

  # Save and return table
  if (!is.null(file)) save_as_docx(resft, path = file)
  
  return(invisible(resft))
}


#--------------------------
# Figure: RR vs PMCI
#--------------------------

#----- Main plot function
plot_RRpred <- function(modlist, var, cities, countrypal, modelpal, ngrid = 100, 
  file = NULL)
{
  
  # Create grid of selected variable
  vargrid <- with(cities, seq(min(cities[,var]), max(cities[,var]), 
    length.out = ngrid))
  
  # Add other variables
  modvars <- lapply(modlist, function(x) all.vars(formula(x$fit))[-1]) |>
    unlist() |> unique() |> setdiff(var)
  othermean <- reframe(cities[, modvars], across(everything(), mean))
  curvedata <- data.frame(vargrid, othermean)
  names(curvedata)[1] <- var
  
  # Predict
  rrpreds <- lapply(modlist, function(x) predict(x$fit, curvedata, ci = T)) |> 
    do.call(what = cbind) |>
    exp() |>
    as.data.frame()
  
  # Names and melt
  allnames <- outer(names(modlist), c("fit", "low", "high"), 
    paste, sep = "_") |> t() |> c()
  colnames(rrpreds) <- allnames
  curvedata <- cbind(curvedata, rrpreds)
  curvedata <- pivot_longer(curvedata, cols = all_of(allnames), 
    names_to = c("model", ".value"), names_sep = "_")

  # Plot
  rrplot <- ggplot(curvedata) + theme_bw() + 
    geom_line(aes(x = .data[[var]], y = fit, col = model), linewidth = 1) + 
    geom_line(aes(x = .data[[var]], y = low, col = model), linetype = 2) + 
    geom_line(aes(x = .data[[var]], y = high, col = model), linetype = 2) +
    scale_color_manual(values = modelpal, breaks = names(modelpal), 
      name = "Model", na.value = "black") +
    geom_hline(yintercept = 1) + 
    labs(x = var, y = "RR")
  
  # If a single model, plot BLUPs
  if (length(modlist) == 1){
    cities[row.names(model.frame(modlist[[1]]$fit)), "blup"] <- 
      exp(blup(modlist[[1]]$fit))
    rrplot <- rrplot + 
      new_scale("color") +
      geom_point(aes(x = .data[[var]], y = blup, col = countryname, 
        shape = countryname), data = cities) +
      scale_color_manual(values = countrypal$color, 
        breaks = names(countrypal$color), 
        name = "Country") + 
      scale_shape_manual(values = countrypal$shape, 
        breaks = names(countrypal$shape), 
        name = "Country")
  }
  
  # Save and return
  if (!is.null(file)) ggsave(file, rrplot, height = 5, width = 8)
  
  return(invisible(rrplot))
}

#----- Function to create palettes (reused elsewhere)
make_palettes <- function(cities){
  
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
  shppal <- rep_len(15:18, sum(!cntrpaldf$isUSA))
  names(shppal) <- subset(cntrpaldf, !isUSA, countryname, drop = T)
  usashp <- rep_len(11:12, length(grep("^USA", cntrpaldf$countryname)))
  names(usashp) <- subset(cntrpaldf, isUSA, countryname, drop = T)
  shppal <- c(shppal, usashp)
  
  # Return
  list(color = cntrpal, shape = shppal)
}
