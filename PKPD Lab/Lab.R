# Clear the environment
rm(list = ls())

# Load required packages
library(lattice)
library(ggplot2)
library(deSolve)

# Set working directory (ensure this path is correct)
setwd("C:/Users/Ethan/Documents/BIOE8510/RSudio/PKPD Lab")

# Load the data
sad <- read.csv("cureaprobex_SAD.csv")

sad <- subset(sad,ID==1)

# Optional: Filter data for a specific subject ID
# sad1 <- subset(sad, ID == 1)

# Plot the raw data
ggplot(sad) +
  geom_path(aes(x = time_hr, y = conc_ng_ml, color = factor(dose), group = ID)) +
  facet_wrap(~ID)

############################################
# Define the two-compartment PK model
############################################

model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dA_depot <- -Ka * A_depot  # Drug leaves the depot
    dA_c <- (Ka * A_depot) - (Ke * A_c) - (K12 * A_c) + (K21 * A_p)  # Central compartment (1)
    dA_p <- (K12 * A_c) - (K21 * A_p)  # Peripheral compartment (2)
    
    #return list of derivatives
    list(c(dA_depot, dA_c, dA_p))  
  })
}

# Set initial conditions and parameters
Dose <- 1e6 

inits <- c(A_depot = Dose, 
           A_c = 0, 
           A_p = 0)

theta <- c(Ka = 0.01, 
           Ke = 9, 
           K12 = 3, 
           K21 = 9)

# Simulation time
tlast <- 24
times <- seq(0, tlast, by = 0.01)

# Run the ODE simulation
X <- ode(y = inits, times = times, func = model, parms = theta)
X <- as.data.frame(X)  # Convert to dataframe for plotting

# Plot both central and peripheral compartments
ggplot() +
  geom_path(data = X, aes(x = time, y = A_c, color = 'Central (A_c)'), size = 1) +
  geom_path(data = X, aes(x = time, y = A_p, color = 'Peripheral (A_p)'), size = 1) +
  geom_point(data = sad, aes(x = time_hr, y = conc_ng_ml, color = 'Observed Data'), size = 2) +
  labs(title = "Two-Compartment PK Model Simulation vs. Observed Data",
       x = "Time (hours)", y = "Concentration (ng/mL)") +
  theme_minimal() +
  scale_color_manual(values = c('Central (A_c)' = 'red', 
                                'Peripheral (A_p)' = 'blue', 
                                'Observed Data' = 'black')) +
  theme(legend.title = element_blank())
