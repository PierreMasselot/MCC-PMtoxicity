#####################################################################
#
#                 MCC HetPoll
#               Part 5: Plots
#
#####################################################################

#----------------------------
# Plot how RR affects CAPI
#----------------------------

#----- Predictions for various values

# Prediction data.frame
capi_pred <- model.frame(capifit) |> 
  mutate(pcs = t(colMeans(pcs)))

# Prediction
capi_RR <- predict(capifit, capi_pred, se = T) |>
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

#----------------------------
# Plot how RR changes when adjusted
#----------------------------

# Plot
ggplot(cities, aes(x = coef, y = coefadj)) + 
  geom_hline(yintercept = 0, col = "grey") + 
  geom_vline(xintercept = 0, col = "grey") + 
  geom_point() + 
  geom_smooth(method = "lm", col = 4) + 
  geom_text(aes(label = cityname), nudge_y = .01,
    data = subset(cities, abs(coef - coefadj) > .1)) + 
  geom_abline(slope = 1, intercept = 0, col = 1) + 
  labs(x = "Unadjusted log(RR)", y = "Adjusted log(RR)") +
  theme_bw()

# Save
ggsave("figures/RRadjustment_scatter.pdf")
