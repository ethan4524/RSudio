#Load packages
library(tidyverse)
library(deSolve)
library(ggpubr)

#===================Case 1: Part E======================#

# Define a range of Fixaprob concentrations (pM)
fixaprob_concentrations <- seq(0, 100, by = 1)  # From 0 to 100 pM
Kd_fixaprob <- 10  # pM

# Calculate fraction bound for each Fixaprob concentration
fb <- fixaprob_concentrations / (Kd_fixaprob + fixaprob_concentrations)

# Create a data frame to store the results
results_case1 <- data.frame(I = fixaprob_concentrations, fb = fb)

# Plot fraction bound vs Fixaprob concentration with a vertical line at Kd = 10 pM
ggplot(results_case1, aes(x = I, y = fb)) +
  geom_line(color = "blue", size = 1.2) +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", size = 1) +
  labs(
    x = "Fixaprob Concentration (pM)",
    y = "Fraction Bound",
    title = "Fraction Bound vs Fixaprob Concentration (Case 1)"
  ) +
  annotate("text", x = 12, y = 0.5, label = "Kd = 10 pM", color = "red", angle = 90, vjust = -0.5) +
  theme_minimal()





