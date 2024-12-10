#Load packages
library(tidyverse)
library(deSolve)
library(ggpubr)

#===================Case 1: Part E======================#

# Define a range of Fixaprob concentrations (pM)
fixaprob_concentrations <- seq(0, 100, by = 1)  # From 0 to 100 pM
Kd <- 10  # pM
Pss=20#pM

Kie_fixaprob = 10

#alculate fraction bound for each Fixaprob conc
fb = Pss / (Pss+Kd*(1+(fixaprob_concentrations / Kie_fixaprob)))
# Create a data frame to store the results
results_case1 <- data.frame(I = fixaprob_concentrations, fb = fb)

# Plot fraction bound vs Fixaprob concentration with a horizontal line at y = 0.222
ggplot(results_case1, aes(x = I, y = fb)) +
  geom_line(color = "blue", size = 1.2) +
  geom_hline(yintercept = 0.222, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 60, color = "red", linetype = "dotted", size = 1) +
  labs(
    x = "Fixaprob Concentration (pM)",
    y = "Fraction Bound",
    title = "Fraction Bound vs Fixaprob Concentration (Case 1)"
  ) +
  annotate("text", x = 12, y =0.225, label = "Fb_min = 0.222", color = "red", vjust = -0.5) +
  theme_minimal()


# Find the index where the fraction bound is closest to 0.222
target_fb <- 0.222
closest_index <- which.min(abs(results_case1$fb - target_fb))

# Get the corresponding Fixaprob concentration
min_conc <- results_case1$I[closest_index]
min_conc

