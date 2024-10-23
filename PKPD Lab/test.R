#Load packages
library(lattice)
library(ggplot2)
library(deSolve)


# Set working directory
setwd("C:/Users/Ethan/Documents/BIOE8510/RSudio/PKPD Lab")

# Read the data
sad = read.csv("cureaprobex_SAD.csv")

ggplot(sad) + geom_path(aes(x=time_hr, y = conc_ng_ml , color = dose, group = ID))
ggplot(sad) + geom_path(aes(x=time_hr, y = conc_ng_ml , color = factor(dose))) + facet_wrap(~ID)

sad=sad
########################################### Define 1 compartment PK Model with oral dosing #####################
# Define the two-compartment model
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Differential equations
    dA_depot <- -Ka * A_depot  # Drug leaves the GI tract
    dA_c <- Ka * A_depot - Ke * A_c - K12 * A_c + K21 * A_p  # Central compartment
    dA_p <- K12 * A_c - K21 * A_p  # Peripheral compartment
    
    # Return the rate of change
    list(c(dA_depot, dA_c, dA_p))
  })
}


# Set initial conditions (amounts in each compartment at t=0)
initial_state <- c(A_depot = 1,  # Initial dose (mg)
                   A_c = 0,        # Central compartment
                   A_p = 0)        # Peripheral compartment

# Set parameters
parameters <- c(Ka = 1.0,    # Absorption rate constant (1/hour)
                Ke = 0.1,     # Elimination rate constant from central (1/hour)
                K12 = 0.05,   # Transfer rate from central to peripheral (1/hour)
                K21 = 0.03)   #    # Transfer rate from peripheral to central (1/hour)

# Define time points for the simulation (0 to 48 hours)
times <- seq(0, 48, by = 0.1)

# Solve the ODEs
out <- ode(y = initial_state, times = times, func = model, parms = parameters)

# Convert the output to a data frame for plotting
out <- as.data.frame(out)

# Plot the results
library(ggplot2)
ggplot(out, aes(x = time)) +
  geom_line(aes(y = A_c, color = "Central Compartment"), size = 1) +
  geom_line(aes(y = A_p, color = "Peripheral Compartment"), size = 1) +
  labs(title = "Two-Compartment Pharmacokinetic Model",
       x = "Time (hours)", y = "Amount of Drug (mg)") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.title = element_blank())