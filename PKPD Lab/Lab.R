# Clear the environment
rm(list = ls())

# Load required packages
library(ggplot2)
library(deSolve)

# Set working directory (modify as needed)
#setwd("C:/Users/Ethan/Documents/BIOE8510/RSudio/PKPD Lab")
setwd("C:/Users/ethan/OneDrive/Documents/RStudio/RSudio/PKPD Lab")

# Load the data
sad <- read.csv("cureaprobex_SAD.csv")

# Filter data for subject ID = 1
sad <- subset(sad, dose == 10)
#sad <- subset(sad, ID == 1)

# Convert concentration ng/ml --> ng/L
sad$conc_ng_L <- sad$conc_ng_ml * 1000

############################################
# Define the two-compartment PK model
############################################

model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    # Depot (GI) compartment dynamics
    dA_depot_dt <- -Ka * a_depot
    dC1_dt <- (Ka * a_depot - k12 * C1 + k21 * C2 - Ke * C1) / V1
    dC2_dt <- k12 * C1 - k21 * C2
    
    # Return list of derivatives
    list(c(dA_depot_dt, dC1_dt, dC2_dt))
  })
}

# Define the initial dose (mg)
Dose <- 1e+7 #10mg in ng

# Initial conditions
inits <- c(a_depot = Dose,  # Amount in the depot
           C1 = 0,           # Initial concentration in central 
           C2 = 0)           # Initial concentration in peripheral 

# Model parameters
theta <- c(Ka = 1,        # (/hr)
           V1 = 86,   # Volume of central compartment (L)
           k12 = 0.05,      #rate from central to peripheral (/hr)
           k21 = 0.1,      #  rate from peripheral to central (/hr)
           Ke = 5)       # (/hr)

# Simulation time 
tlast <- 48
times <- seq(0, tlast, by = 0.01)

# Run the ODE simulation
X <- ode(y = inits, times = times, func = model, parms = theta)
X <- as.data.frame(X)  # Convert to dataframe for plotting

# Plot peripheral and central concentrations with observed data
ggplot() +
  # Observed data points
  geom_point(data = sad, aes(x = time_hr, y = conc_ng_L, color = "Observed"), size = 2) +
  # Central compartment concentrations
  geom_line(data = X, aes(x = time, y = C1, color = "Central Compartment (C1)"), size = 1) +
  # Peripheral compartment concentrations
  geom_line(data = X, aes(x = time, y = C2, color = "Peripheral Compartment (C2)"), size = 1) +
  labs(title = "Peripheral & Central Compartment vs Observed Concentration",
       x = "Time (hr)", y = "Concentration (ng/L)", color = "Legend") +
  theme_minimal()




################### Optimization using Optim ########################

#Define an objective function. This function should take in a list of parameters as input and return the
#objective function value

objfxn <- function(theta) {
  inits <- c(a_depot = Dose, C1 = 0, C2 = 0)
  times <- sad$time_hr  # Use observed time points for simulation
  #times <- rep(times,8)
  
  # Run the model simulation
  out <- ode(y = inits, times = times, func = model, parms = theta)
  out<-rep(out,8)
  
  out <- as.data.frame(out)  # Convert to data frame for easier manipulation
  # Extract predicted concentrations from the central compartment (C1)
  pred_c1 <- out$C1
  
  # Calculate the mean squared error (MSE) between observed and predicted concentrations
  mse <- mean((sad$conc_ng_L - pred_c1)^2)
  
  # Print for debugging (optional)
  print(mse)
  print(theta)
  
  # Return the MSE as the objective function value
  return(mse)
}


#test objfxn

#Test your objective function by running with a test case. It should return a single value - the MSE.
objfxn(theta) #Running this line should return a single value for your objective function.



fit = optim(theta, objfxn, method = "BFGS", hessian = T)


#Results of the optimization are stored in the variable "fit"
fit



#objm time= x- sse = 
# use the rep functin to repeat the times 0->time_end for each subject
