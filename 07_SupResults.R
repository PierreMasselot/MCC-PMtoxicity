#####################################################################
#
#                 MCC HetPoll
#               Part 7: Supplementary plots
#
#####################################################################

#----------------------------
# Residuals
#----------------------------

## Reuses palettes from Figure 1

# Plot association
residplot <- ggplot(cities) + theme_bw() + 
  geom_point(aes(x = PMCI, y = exp(residuals), col = countryname, 
    shape = countryname)) + 
  geom_hline(yintercept = 1) + 
  scale_color_manual(values = cntrpal, breaks = names(cntrpal), 
    name = "Country") + 
  scale_shape_manual(values = shppal, breaks = names(cntrpal), 
    name = "Country") + 
  labs(x = "PMCI", y = "log(RR) residual")

# Save
ggsave("figures/SupFig1_Residuals.png")
