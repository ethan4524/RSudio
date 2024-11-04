# Load necessary libraries
library(deSolve)
library(ggplot2)

# Load pharmacodynamic data
mad_PD <- read.csv("cureaprobex_MAD_PD.csv")
Dose = 1
mad_PD <- subset(mad_PD, dose == Dose)

# Define parameters (you might want to adjust these)
theta <- list(
  Ka = 0.1,          # absorption rate constant
  V1 = 10,          # volume of the central compartment
  k12 = 0.5,        # rate constant for distribution from central to peripheral
  k21 = 0.3,        # rate constant for distribution from peripheral to central
  Ke = 0.05,        # elimination rate constant
  k_inhibition = 0.01,  # rate of inhibition of Wrkin1
  IC50 = 5,         # IC50 for inhibition
  Imax = 1)

# Define initial conditions
Dose <- 1  # 1 mg = 1000 µg
inits <- c(a_depot = Dose, 
           C1 = 0, 
           C2 = 0,
           Wrkin1 = 3.106001)  # Initial Wrkin1 concentration based on observed data

# Simulation settings
hours <- 24          # Simulate for 24 hours
times <- seq(0, hours, by = 1)  # Time vector in hours

# Define the model function
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    # Depot (GI) compartment dynamics
    dA_depot_dt <- -Ka * a_depot
    dC1_dt <- (Ka * a_depot - k12 * C1 + k21 * C2 - Ke * C1) / V1
    dC2_dt <- k12 * C1 - k21 * C2
    
    # Inhibition of Wrkin1 production
    inhibition <- 1 - ((Imax * C1) / (IC50 + C1))
    
    # Rate of change of Wrkin1
    dWrkin1_dt <- -k_inhibition * Wrkin1 * inhibition
    
    # Return list of derivatives
    list(c(dA_depot_dt, dC1_dt, dC2_dt, dWrkin1_dt))
  })
}

# Run the simulation using ode solver
out <- ode(y = inits, times = times, func = model, parms = theta)

# Convert output to a data frame for plotting
results_df <- as.data.frame(out)

# Plot observed data and simulated response
ggplot() +
  geom_line(data = mad_PD, aes(x = time_hr, y = wrkin1_pg_ml), color = "blue", size = 1, alpha = 0.7) +  # Observed data
  geom_line(data = results_df, aes(x = time, y = Wrkin1), color = "red", size = 1, linetype = "dashed") +  # Simulated response
  labs(title = "Observed vs. Simulated Wrkin1 Concentration",
       x = "Time (hours)",
       y = "Wrkin1 Concentration (pg/mL)",
       caption = "Blue: Observed Data; Red: Simulated Response") +
  theme_minimal() +
  xlim(0, 24)  # Limit x-axis to first 24 


# Example: Adjust parameters and run simulations
adjusted_params <- function(ka, k_inh, ic50, imax) {
  theta <- list(
    Ka = ka,
    V1 = 10,
    k12 = 0.5,
    k21 = 0.3,
    Ke = 0.05,
    k_inhibition = k_inh,
    IC50 = ic50,
    Imax = imax
  )
  
  out <- ode(y = inits, times = times, func = model, parms = theta)
  return(as.data.frame(out))
}

# Example usage
results_df_adjusted <- adjusted_params(0.15, 0.02, 4, 1)  # Adjust these values as needed

# Plot adjusted results
ggplot() +
  geom_line(data = mad_PD, aes(x = time_hr, y = wrkin1_pg_ml), color = "blue", size = 1, alpha = 0.7) +  # Observed data
  geom_line(data = results_df_adjusted, aes(x = time, y = Wrkin1), color = "red", size = 1, linetype = "dashed") +  # Adjusted simulated response
  labs(title = "Adjusted Simulated Wrkin1 Concentration",
       x = "Time (hours)",
       y = "Wrkin1 Concentration (pg/mL)",
       caption = "Blue: Observed Data; Red: Adjusted Simulated Response") +
  theme_minimal() +
  xlim(0, 24)  # Limit x-axis to first 24 hours
