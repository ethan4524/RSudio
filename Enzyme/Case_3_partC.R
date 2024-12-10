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
    PR_ss = (Etot_ss * Pss) / (Kd + Pss)
    
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
        Iblocextra = 15,      # Initial Blocextra concentration
        Pss = Pss,
        n = 3,               # Feedback exponent
        Kd = Kd
      )

# Time parameters
tlast = 60 * 60 * 24 * 3  # 3 days in seconds
times <- seq(0, tlast, 60)  # Every 60 seconds

# Solve the differential equations using ode function
X <- data.frame(ode(inits, times, model, theta))

# Calculate the fraction bound
X$fb = X$P / (theta["Kd"] + X$P)

# Reshape the data into long format for easy plotting with ggplot
X_long <- X %>%
  pivot_longer(cols = c(Et, St, P, fb), names_to = "Variable", values_to = "Value")

# Create a combined plot with 4 subplots and legend
ggplot(X_long, aes(x = time, y = Value, color = Variable)) +
  geom_line(size = 1) +
  labs(x = "Time (seconds)", y = "Concentration / Fraction Bound", title = "Enzyme Kinetics Over Time") +
  theme_minimal() +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +  # Create separate subplots for each variable
  scale_color_manual(values = c("blue", "red", "green", "purple"), 
                     labels = c("Enzyme (Et)", "Substrate (St)", "Product (P)", "Fraction Bound")) +  # Customize colors and legend labels
  theme(legend.title = element_blank(), 
        legend.position = "bottom")

# Get the final concentration of free Probletide (Et)
final_Et <- tail(X$Et, 1)

# Print the final concentration of Et
final_Et