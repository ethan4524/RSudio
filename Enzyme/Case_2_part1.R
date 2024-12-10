# Load packages
rm(list=ls())

library(tidyverse)
library(deSolve)
library(ggpubr)

# Specify deSolve Model
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    # Enzyme-substrate complex with Fixaprob and Blocextra as competitive inhibitors
    ES = Et * St / ((Km * (1 + (Iblocextra / Kie_blocextra))) + St)
    
    # Free enzyme
    E = Et - ES
    
    # Free substrate
    #S = St - ES
    
    # Fraction bound
    fb = P / (P+Kd*(1+(Ifixaprob / Kie_fixaprob)))
    
    # Differential equations
    dEt = (Prode *  (fbss/fb)^n)- Kdege*E
    dSt = Prods - Kdegs * St - Kcat * ES
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
  Ifixaprob = 60,       # Initial Fixaprob concentration
  Iblocextra = 0,      # Initial Blocextra concentration
  Pss = Pss,
  n = 4.4,               # Feedback exponent
  Kd = Kd,
  fbss = fbss
)

# Simulation time (3 days in seconds)
tlast = 60 * 60 * 24 * 3  # 3 days in seconds
times <- seq(0, tlast, 60)  # Every 60 seconds

# Simulate the model
X <- data.frame(ode(inits, times, model, theta))

# Calculate fraction bound (use time column from X)
X <- X %>%
  mutate(fb = P / (Kd + P))  # Use Kd from theta

# Convert data to long format for ggplot
X_long <- X %>%
  pivot_longer(cols = c(Et, St, P, fb), 
               names_to = "Variable", 
               values_to = "Value")

# Plot with facet_wrap
ggplot(X_long, aes(x = time, y = Value)) +
  geom_line(color = "steelblue", size = 1) +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(x = "Time (seconds)", y = "Concentration", title = "Simulation of Et, St, P, and Fraction Bound") +
  theme_minimal()



#=======================Part A================================

# Function to calculate fraction bound and enzyme concentration
calculate_fraction_bound <- function(fixaprob_concentration) {
  # Update Ifixaprob value in the theta list
  theta["Ifixaprob"] <- fixaprob_concentration
  # Simulate the model
  X <- data.frame(ode(inits, times, model, theta))
  
  # Calculate fraction bound
  X$fb <- X$P / (X$P+theta["Kd"]*(1+(theta["Ifixaprob"] / theta["Kie_fixaprob"])))
  
  # Get final fraction bound value (last time point)
  final_fb <- tail(X$fb, n = 1)
  
  # Get final Enzyme concentration value (last time point)
  final_Et <- tail(X$Et, n = 1)
  print(X$Et)
  # Return as a data frame row
  print(data.frame(fixaprob_concentration, final_fb, final_Et))
  return(data.frame(fixaprob_concentration, final_fb, final_Et))
}

# Create a sequence of Fixaprob concentrations to test
fixaprob_concentrations <- seq(0, 1000, by = 5)

#runs the function above with the list of fixaprob concentrations as inputs
results_list <- lapply(fixaprob_concentrations, calculate_fraction_bound)
#store results in dataframe
results <- bind_rows(results_list)

# Print the results
print(results)

# Plot
ggplot(results, aes(x = final_fb, y = final_Et)) +
  geom_line(color = "purple", size = 1) +
  geom_point(color = "purple", size = 2) +
  labs(
    x = "Final Fraction Bound (Fb)",
    y = "Final Enzyme Concentration (Et)",
    title = "Enzyme Concentration vs Fraction Bound"
  ) +
  theme_minimal()
