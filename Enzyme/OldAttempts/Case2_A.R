#Load packages
library(tidyverse)
library(deSolve)
library(ggpubr)

# ------------------------ LEGEND FOR VARIABLES ------------------------

# Et: Enzyme concentration (Extrase)
#      - Represents the total concentration of the enzyme involved in the reaction.

# St: Substrate concentration (Substrase)
#      - Represents the total concentration of the substrate available for the enzyme to act upon.

# ES: Enzyme-substrate complex
#      - Calculated internally in the model as the complex formed between the enzyme (Extrase) and the substrate (Substrase).

# E: Free enzyme
#      - The enzyme that is not bound to the substrate.

# P: Product concentration (Problatide)
#      - Represents the total concentration of the product formed by the reaction.

# I: Inhibitor concentration (Fixaprob)
#      - Represents the concentration of the inhibitor used to block enzyme activity.

# Km: Michaelis constant
#      - A measure of the substrate concentration at which the reaction rate is half-maximal.

# Kcat: Catalytic rate constant
#      - The number of substrate molecules converted to product per enzyme molecule per second.

# Kdege: Degradation rate constant for the enzyme
#      - The rate at which the enzyme degrades over time.

# Kdegs: Degradation rate constant for the substrate
#      - The rate at which the substrate degrades over time.

# Prode: Enzyme production rate
#      - The rate at which the enzyme is produced.

# Prods: Substrate production rate
#      - The rate at which the substrate is produced.

# Kdegp: Degradation rate constant for the product
#      - The rate at which the product degrades over time.

# ------------------------ END OF LEGEND ------------------------

model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {  # Get named contents of inits and theta
    
    # Calculate enzyme-substrate complex and free enzyme
    ES = Et * St / (Km * (1 + Ie / Kie) + St)
    E = Et - ES
    S = St - ES
    
    # Fraction bound with inhibitor (Fixaprob)
    fb = P / (Kd + P) * (1 / (1 + Fixaprob / Kfix))
    
    # Differential equations
    dEt = (Prode * (fbss / fb)^n) - Kdege * E
    dSt = Prods - Kdegs * S - Kcat * ES
    dP = Kcat * ES - Kdegp * P
    
    # Return list of derivatives
    list(c(dEt, dSt, dP))
  })
}

###################################################3
# Known constants
Km = 1.25 * 1000 # pM
Kcat = 4 # /sec
Kdege = 1 # /sec (Ke)
Kdegs = 0.1 # /sec (the substrate can only be eliminated by converting into product)
Kd = 10 # pM
Kie = 9500 # pM (Inhibitor binding affinity for the enzyme)

# Concentrations measured at steady state
Pss = 20 # pmol/L
Sss = 10 # pmol/L
Etot_ss = 0.5 # pmol/L

# Calculated concentrations at SS

# Enzyme-substrate complex at steady-state
ES_ss = Etot_ss * Sss / (Km + Sss)

# Free enzyme at steady state
E_ss = Etot_ss - ES_ss

fbss = Pss / (Pss + Kd)

# Calculated production and elimination rates 
Prode = Kdege * E_ss    # Must equal free-enzyme elimination rate
Prods = Kdegs * Sss + Kcat * ES_ss  # Must equal free enzyme elimination rate + product formation rate
Kdegp = Kcat * ES_ss / Pss  # Calculate from steady-state condition for Product

#####################################3
inits = c(Et = Etot_ss,
          St = Sss + ES_ss,
          P = Pss
          )  # Initialize Fixaprob as 0 for Case 2

theta = c(Km = Km,
          Kcat = Kcat,
          Kdege = Kdege,
          Kdegs = Kdegs,
          Prode = Prode,
          Prods = Prods,
          Kdegp = Kdegp,
          n = 3,
          Kd = Kd,
          Ie = 100,  # Initial concentration of the inhibitor (Fixaprob)
          Kie = Kie,
          Fixaprob = 0,  # Initial concentration of Fixaprob (change as needed)
          Kfix = 1)  # Binding affinity for Fixaprob

# Simulation time
tlast = 60 * 24 * 3  # 3 days in minutes
times <- seq(0, tlast, 60)  # Simulate every minute

# Simulate
out <- ode(inits, times, model, theta)

# Convert to data frame and include fb
X <- as.data.frame(out)

# Calculate fraction bound
X$fb = X$P / (theta["Kd"]  + X$P)

# Post-hoc calculations

# Reshape data to long format for faceting
X_long <- X %>%
  pivot_longer(cols = c(Et, St, P, fb), 
               names_to = "Variable", 
               values_to = "Value")

# Plot using facet_wrap
ggplot(X_long, aes(x = time, y = Value)) +
  geom_line(color = "blue") +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "Simulation of Et, St, P, and Fraction Bound",
       x = "Time (seconds)",
       y = "Concentration or Fraction Bound") +
  theme_minimal()
