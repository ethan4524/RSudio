# Load necessary packages
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
Kdegs = 0.1  # /sec
Kd = 10  # pM

# Concentrations measured at steady state
Pss = 20  # pmol/L
Sss = 10  # pmol/L
Etot_ss = 0.5  # pmol/L

# Calculate steady-state concentrations
ES_ss = Etot_ss * Sss / (Km + Sss)
E_ss = Etot_ss - ES_ss
Prode = Kdege * E_ss
Prods = Kdegs * Sss + Kcat * ES_ss
Kdegp = Kcat * ES_ss / Pss
PR_ss = (Etot_ss * Pss) / (Kd + Pss)

# Initial conditions for simulation
inits = c(Et = Etot_ss, St = Sss + ES_ss, P = Pss)

# Parameter list
theta <- c(
  Km = Km,
  Kcat = Kcat,
  Kdege = Kdege,
  Kdegs = Kdegs,
  Prode = Prode,
  Prods = Prods,
  Kdegp = Kdegp,
  Kie_fixaprob = 10,  # pM
  Kie_blocextra = 1,  # pM
  Ifixaprob = 0,  # Initial Fixaprob concentration (to be varied)
  Iblocextra = 1953,  # Initial Blocextra concentration
  Pss = Pss,
  n = 3,  # Feedback exponent
  Kd = Kd
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
  labs(x = "Time (minutes)", y = "Concentration", title = "Simulation of Et, St, P, and Fraction Bound") +
  theme_minimal()



