# Load packages
library(tidyverse)
library(deSolve)
library(ggpubr)

# Clear the workspace
rm(list = ls())

# Specify the deSolve Model
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    # Enzyme-substrate complex with Blocextra (enzyme inhibitor)
    ES = Et * St / (Km * (1 + Ib / Kb) + St)
    E = Et - ES
    S = St - ES
    
    # Fraction bound with Fixaprob (receptor inhibitor)
    fb = P / (Kd + P) * (1 / (1 + Fixaprob / Kfix))  # Adjusting feedback for receptor inhibition
    
    # Differential equations
    dEt = (Prode * (fbss / fb)^n) - Kdege * E
    dSt = Prods - Kdegs * S - Kcat * ES
    dP = Kcat * ES - Kdegp * P
    
    # Return list of derivatives
    list(c(dEt, dSt, dP))
  })
}

# Known constants
Km <- 1.25 * 1000    # nM 
Kcat <- 4            # /sec
Kdege <- 1           # /sec
Kdegs <- 0.1         # /sec (the substrate can only be eliminated by converting into product)
Kd <- 10              # pM (binding affinity for feedback)

# Concentrations measured at steady state
Pss <- 20            # pM
Sss <- 10            # pM
Etot_ss <- 0.5       # pM

# Calculated concentrations at steady state
ES_ss <- Etot_ss * Sss / (Km + Sss)  # Enzyme-substrate complex at steady state
E_ss <- Etot_ss - ES_ss              # Free enzyme at steady state

# Calculated production and elimination rates
Prode <- Kdege * E_ss                      # Must equal free-enzyme elimination rate
Prods <- Kdegs * Sss + Kcat * ES_ss        # Must equal substrate degradation + product formation rate
Kdegp <- Kcat * ES_ss / Pss                # Product degradation rate

# Fraction bound at steady state
fbss <- Pss / (Pss + Kd)

# Initial conditions
inits <- c(
  Et = Etot_ss,
  St = Sss + ES_ss,
  P = Pss
)

# Parameters
theta <- c(
  Km = Km,             # Michaelis constant
  Kcat = Kcat,         # Catalytic rate constant
  Kdege = Kdege,       # Enzyme degradation rate constant
  Kdegs = Kdegs,       # Substrate degradation rate constant
  Prode = Prode,       # Enzyme production rate
  Prods = Prods,       # Substrate production rate
  Kdegp = Kdegp,       # Product degradation rate constant
  n = 3,               # Feedback exponent
  Kd = Kd,             # Dissociation constant for feedback
  Fixaprob = 0,        # Fixaprob concentration (set to 0 for Case 3)
  Kfix = 0,            # Placeholder binding affinity for Fixaprob (not used here)
  Ib = 0,              # Initial Blocextra concentration (to be varied in simulation)
  Kb = 1               # Binding affinity of Blocextra (1 pM as per instructions)
)

# Simulation time: 3 days (in seconds)
tlast <- 60 * 60 * 24 * 3  # 3 days
times <- seq(0, tlast, by = 60)  # Every 60 seconds (1 minute)

# Define a range of Blocextra concentrations
blocextra_concentrations <- seq(0, 100, by = 5)  # From 0 to 100 pM in steps of 5

# Create an empty data frame to store results
results <- data.frame()

# Loop over each concentration of Blocextra
for (Ib_value in blocextra_concentrations) {
  
  # Update Blocextra concentration in theta
  theta["Ib"] <- Ib_value
  
  # Simulate the model
  out <- ode(inits, times, model, theta)
  
  # Convert output to data frame and add Blocextra concentration
  X <- as.data.frame(out)
  X$Ib <- Ib_value
  
  # Calculate fraction bound
  X$fb <- X$P / (theta["Kd"] + X$P)
  
  # Append results
  results <- rbind(results, X)
}

# Average the fraction bound for each Blocextra concentration
avg_results <- results %>%
  group_by(Ib) %>%
  summarize(mean_fb = mean(fb, na.rm = TRUE))

# Plot Fraction Bound vs Blocextra Concentration
ggplot(avg_results, aes(x = Ib, y = mean_fb)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(title = "Fraction Bound vs Blocextra Concentration",
       x = "Blocextra Concentration (pM)",
       y = "Mean Fraction Bound (fb)") +
  theme_minimal()
