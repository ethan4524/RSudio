
#Part 6:
##========#========#========#========#========#========#========#========
#========Simulating the Pharmacodynamic response.#==========#============
##========#========#========#========#========#========#========#========

# Step 1: Observe the data
# Load data
mad_PD <- read.csv("cureaprobex_MAD_PD.csv")
head(mad_PD)  
Dose <- 1
mad_PD <- subset(mad_PD, dose == Dose)  

# Filter observed data to only include the first 24 hours
#mad_PD_24hr <- subset(mad_PD, time_hr <= 24)

# Define the model with Emax effect on Wrkin1
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    # Depot (GI) compartment dynamics
    dA_depot_dt <- -Ka * a_depot
    dC1_dt <- (Ka * a_depot - k12 * C1 + k21 * C2 - Ke * C1) / V1
    dC2_dt <- k12 * C1 - k21 * C2
    
    # Emax model for pharmacodynamic effect on Wrkin1
    Effect <- Emax * C1 / (EC50 + C1)  # This reflects the effect on production, no additive 1 here
    
    # Rate of change of Wrkin1 based on the pharmacodynamic effect
    dWrkin1_dt <- -k_inhibition * Wrkin1 * Effect
    
    # Return list of derivatives
    list(c(dA_depot_dt, dC1_dt, dC2_dt, dWrkin1_dt))
  })
}

# Define parameters for the model using Emax
theta <- list(
  Ka = fit$par[1],          # absorption rate constant
  V1 = fit$par[2],          # volume of the central compartment
  k12 = fit$par[3],         # rate constant for distribution from central to peripheral
  k21 = fit$par[4],         # rate constant for distribution from peripheral to central
  Ke = fit$par[5],          # elimination rate constant
  Emax = 0.3,               # maximum effect for inhibition
  EC50 = 0.015,               # concentration at which half-maximal effect is achieved
  k_inhibition = 0.1        # rate of inhibition of Wrkin1
)

# Define initial conditions
inits <- c(a_depot = Dose, 
           C1 = 0, 
           C2 = 0,
           Wrkin1 = 3)  # Set initial Wrkin1 concentration based on the observed data

# Simulation settings
days <- 7          # Simulate for 24 hours
times <- seq(0, 24 * days, by = 1)  # Time vector in hours

# Run the simulation using ode solver
out <- ode(y = inits, times = times, func = model, parms = theta)

# Convert output to a data frame
results_df <- as.data.frame(out)

# Plot the Effect over time
ggplot(results_df, aes(x = time, y = Wrkin1)) +
  geom_line(color = "purple", size = 1) +
  labs(title = "Simulated Wrkin1 Concentration Over Time",
       x = "Time (hours)",
       y = "Wrkin1 Concentration (pg/mL)") +
  theme_minimal()

# Overlay observed data for comparison
ggplot() +
  # Observed data for the first 24 hours
  geom_line(data = mad_PD_24hr, aes(x = time_hr, y = wrkin1_pg_ml), color = "blue", size = 1, alpha = 0.7) +
  geom_point(data = mad_PD_24hr, aes(x = time_hr, y = wrkin1_pg_ml), color = "blue", size = 2, alpha = 0.7) +
  
  # Simulated Wrkin1 response for 24 hours
  geom_line(data = results_df, aes(x = time, y = Wrkin1), color = "red", size = 1, linetype = "dashed") +
  
  # Labels and titles
  labs(title = "Observed vs. Simulated Wrkin1 Concentration Over Time (First 24 Hours)",
       x = "Time (hours)",
       y = "Wrkin1 Concentration (pg/mL)",
       caption = "Blue: Observed Data; Red: Simulated Response (Dashed Line)") +
  
  # Theme and customization
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +  # Set x-axis breaks for 24 hours
  xlim(0, days * 24)  # Limit x-axis to the first 24 hours