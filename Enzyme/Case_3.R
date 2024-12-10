# Load packages
rm(list=ls())

library(tidyverse)
library(deSolve)
library(ggpubr)

# Specify deSolve Model
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    # Enzyme-substrate complex with Fixaprob and Blocextra as competitive inhibitors
    ES = Et * St / ((Km * (1 + (Ifixaprob / Kie_fixaprob) + (Iblocextra / Kie_blocextra))) + St)
    
    # Free enzyme
    E = Et - ES
    
    # Free substrate
    S = St - ES
    
    # Fraction bound
    fb = P / (Kd + P)
    
    # Product-receptor complex
    PR = Et * P / (Kd + P)
    
    # Feedback term
    feedback = (PR_ss / PR)^n
    
    # Differential equations
    dEt = (Prode * feedback) - Kdege * E
    dSt = Prods - Kdegs * S - Kcat * ES
    dP = Kcat * ES - Kdegp * P
    
    # Return list of derivatives
    list(c(dEt, dSt, dP))
  })
}
# Known constants
Km = 1.25 * 1000  # pM
Kcat = 4  # /sec
Kdege = 1  # /sec
Kdegs = 0.1  # /sec (the substrate can only be eliminated by converting into product)
Kd = 10  # pM

# Concentrations measured at steady state
Pss = 20  # pmol/L
Sss = 10  # pmol/L
Etot_ss = 0.5  # pmol/L

# Calculated concentrations at SS

# Enzyme substrate complex at steady-state
ES_ss = Etot_ss * Sss / (Km + Sss)

# Free enzyme at steady state
E_ss = Etot_ss - ES_ss

# Calculated production and elimination rates
Prode = Kdege * E_ss  # Must equal free-enzyme elimination rate
Prods = Kdegs * Sss + Kcat * ES_ss  # Must equal free enzyme elimination rate + product formation rate
Kdegp = Kcat * ES_ss / Pss  # Calculate from steady-state condition for Product

fbss = Pss / (Pss + Kd)

PR_ss = (Etot_ss * Pss) / (Kd + Pss)

# Initial conditions
inits = c(Et = Etot_ss,
          St = Sss + ES_ss,
          P = Pss)

# Parameters (theta)
theta <- c(
  Km = Km,
  Kcat = Kcat,
  Kdege = Kdege,
  Kdegs = Kdegs,
  Prode = Prode,
  Prods = Prods,
  Kdegp = Kdegp,
  Kie_fixaprob = 10,  # pM (Fixaprob binding affinity)
  Kie_blocextra = 1,  # pM (Blocextra binding affinity)
  Ifixaprob = 0,       # Initial Fixaprob concentration
  Iblocextra = 0,      # Initial Blocextra concentration
  Pss = Pss,
  n = 3,               # Feedback exponent
  Kd = Kd
)

# Time parameters
tlast = 60 * 60 * 24 * 3  # 3 days in seconds
times <- seq(0, tlast, 60)  # Every 60 seconds

# Define a function to calculate fraction bound for various blocextra concentrations
calculate_fraction_bound <- function(fixaprob_concentration) {
  # Update Iblocextra value in the theta list
  theta["Iblocextra"] <- fixaprob_concentration
  
  # Simulate the model
  X <- data.frame(ode(inits, times, model, theta))
  
  # Calculate fraction bound
  X$fb = X$P / (theta["Kd"] + X$P)
  
  # Get final fraction bound value (last time point)
  final_fb <- tail(X$fb, n = 1)
  
  return(final_fb)
}

# Create a sequence of Blocextra concentrations to test
blocextra_concentrations <- seq(0, 3000, by = 100)  # From 0 to 20 pM

# Calculate fraction bound for each concentration
fraction_bound_values <- sapply(blocextra_concentrations, calculate_fraction_bound)

# Plot fraction bound as a function of Blocextra concentration
result_df <- data.frame(Blocextra = blocextra_concentrations, Fraction_Bound = fraction_bound_values)

ggplot(result_df, aes(x = Blocextra, y = Fraction_Bound)) +
  # Add light green shaded region between x = 2.84 and x = 75 (behind the data)
  geom_rect(aes(xmin = 2.84, xmax = 75, ymin = -Inf, ymax = Inf), 
            fill = "lightgreen", alpha = 0.0) +  # Set the transparency with alpha
  geom_line(color = "steelblue", size = 1) +   # Line for fraction bound
  geom_point(color = "red") +                   # Points for each Blocextra concentration
  geom_hline(yintercept = 0.222, linetype = "dashed", color = "black") +  # Horizontal line at y = 0.222
  labs(x = "Blocextra Concentration (pM)", 
       y = "Fraction Bound", 
       title = "Fraction Bound vs. Blocextra Concentration") +
  theme_minimal()

# Find the index where the fraction bound is closest to 0.222
target_fb <- 0.222
closest_index <- which.min(abs(result_df$Fraction_Bound - target_fb))

# Get the corresponding Blocextra concentration
closest_blocextra_concentration <- result_df$Blocextra[closest_index]
closest_blocextra_concentration
