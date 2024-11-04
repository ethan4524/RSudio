# Clear the environment
rm(list = ls())

# Load required packages
library(ggplot2)
library(deSolve)
library(tidyverse)

# Set working directory (modify as needed)
setwd("C:/Users/Ethan/Documents/BIOE8510/RSudio/PKPD Lab")
#setwd("C:/Users/ethan/OneDrive/Documents/RStudio/RSudio/PKPD Lab")

#Part 1 & 2:
##========#========#========#========#========#========#========#========
#======== Modeling Oral Dosing and  PK Parameter Estimation.#==========
#========#========#========#========#========#========#========#=========


# Load the data
sad <- read.csv("cureaprobex_SAD.csv")

# Filter data for subject ID = 1
sad <- subset(sad, dose == 1)
#sad <- subset(sad, ID == 1)

# Convert concentration ng/ml --> ng/L
sad$conc_ng_L <- sad$conc_ng_ml * 1000
 
############################################
# Define the two-compartment PK model for oral dosing
############################################

model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    dA_depot_dt <- -Ka * a_depot
    dC1_dt <- (Ka * a_depot - k12 * C1 + k21 * C2 - Ke * C1) / V1
    dC2_dt <- k12 * C1 - k21 * C2
    
    # Return list of derivatives
    list(c(dA_depot_dt, dC1_dt, dC2_dt))
  })
}

# Define the initial dose (mg)
Dose <- 1e+7 #10mg in ng
Dose <- Dose / 10

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
  geom_point(data = sad, aes(x = time_hr, y = conc_ng_L, color = "Observed"), size = 2) +
  
  # Central compartment concentrations
  geom_line(data = X, aes(x = time, y = C1, color = "Central Compartment (C1)"), size = 1) +
  
  # Peripheral compartment concentrations
  geom_line(data = X, aes(x = time, y = C2, color = "Peripheral Compartment (C2)"), size = 1) +
  labs(title = "Peripheral & Central Compartment vs Observed Concentration",
       x = "Time (hr)", y = "Concentration (ng/L)", color = "Legend") +
  theme_minimal()


#fixing the times:
unique_time <- unique(sad$time_hr)


#Part 3:
##========#========#========#========#========#========#========#========
#======== Optimization using Optim#=======#========#========#===========
#========#========#========#========#========#========#========#=========

#Define an objective function. This function should take in a list of parameters as input and return the
#objective function value

objfxn <- function(theta) {
  inits <- c(a_depot = Dose, C1 = 0, C2 = 0)
  
  times <- unique_time
  #times <- rep(times,8)
  
  # Run the model simulation
  out <- ode(y = inits, times = times, func = model, parms = theta)
  #out<-rep(out,8) #I'm not sure if this is needed anymore
  
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

#================

# Extracting estimates and standard errors
estimates <- fit$par  # Estimated parameters
se <- sqrt(diag(solve(fit$hessian)))

# Calculate percent standard error
percent_se <- (se / estimates) * 100

# Create a summary table
results_table <- data.frame(
  Parameter = names(estimates),
  Estimate = estimates,
  SE = se,
  Percent_SE = percent_se
)

# Print the results table
print(results_table)

#==================


#Update theta with your optimized parameters
fit$par
theta$Ka = fit$par[1]
theta$V1 = fit$par[2]
theta$k12 = fit$par[3]
theta$k21 = fit$par[4]
theta$Ke = fit$par[5]


#Check 1: did it converge?
fit$convergence   #If this is anything but zero, the optmization did not converge - it was not able to find the minimum value for the objective function

#Check 2: Visual Check. Plot observed and predicted data
theta <- fit$par

# Run the model simulation with optimized parameters
out <- ode(y = inits, times = unique_time, func = model, parms = theta)
out <- as.data.frame(out)

#Calculate model output with final parameter values
sad$Rmod_fit <- approx(out$time, out$C1, xout = sad$time_hr)$y

# Plot observed and predicted values over time
ggplot() +
  geom_point(data = sad, aes(x = time_hr, y = conc_ng_L, color = "Observed")) +
  geom_line(data = sad, aes(x = time_hr, y = Rmod_fit, color = "Predicted")) +
  labs(title = "Observed vs Predicted Concentration Over Time",
       x = "Time (hr)", y = "Concentration (ng/L)", color = "Legend") +
  theme_minimal()



  
  
#Check 3: Residuals - The residuals are the differences between observed and predicted values. 
sad$resid <- sad$conc_ng_L - sad$Rmod_fit

# Histogram of residuals
hist(sad$resid, main = "Histogram of Residuals", xlab = "Residuals")

# Residuals vs. Time Plot a horizontal reference line at 0
ggplot(sad, aes(x = time_hr, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs. Time", x = "Time (hr)", y = "Residuals") +
  theme_minimal()

#Check 4: Observed vs. Predicted
# - add a diagonal reference line with a slope of 1
ggplot(sad, aes(x = Rmod_fit, y = conc_ng_L)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Observed vs. Predicted Concentrations",
       x = "Predicted Concentration (ng/L)", y = "Observed Concentration (ng/L)") +
  theme_minimal()

#notes:
#Rmax --> max R as t goes to infinity
#k --> growth rate constant
#t50 --> time of half maximum response




#Part 5:
##========#========#========#========#========#========#========#========
#========Evaluate predictive ability for multiple oral dosing.#==========
#========#========#========#========#========#========#========#=========

############################################
# Define the two-compartment PK model again for my sanity
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

#update theta (just to be sure)
theta$Ka = fit$par[1]
theta$V1 = fit$par[2]
theta$k12 = fit$par[3]
theta$k21 = fit$par[4]
theta$Ke = fit$par[5]  

#create list of doses to loop through 
Doses <- c(0.1,1,10) #mg

#nitialize an empty data frame to store results for ALL doses
all_results <- data.frame()

#loop through each dose
for (Dose in Doses) {
  
  # Read in dataset and filter for the current dose
  mad <- read.csv("cureaprobex_MAD_PK.csv")
  mad <- subset(mad, dose == Dose)
  
  # Convert concentration ng/ml to ng/L (I think this is correct)
  mad$conc_ng_L <- mad$conc_ng_ml * 1000
  #mad$conc_ng_L <- mad$conc_ng_ml
  
  #time stuff for simulation
  days <- 24
  times <- seq(0, 24, by = 1/60)  
  
  #initial conditions
  inits <- c(a_depot = Dose,
             C1 = 0, 
             C2 = 0)  
  dose_results <- data.frame()  # To store results for the current dose
  
  # Simulate each dose per day, store data (this is similar to the code provided)
  for (i in 1:days) {
    X <- data.frame(ode(inits, times, model, theta))
    
    # Update initial conditions for the next day
    inits["C1"] <- tail(X$C1, n = 1)
    inits["C2"] <- tail(X$C2, n = 1)
    inits["a_depot"] <- tail(X$a_depot, n = 1) + Dose  # Add new dose each day
    
    # Update the time to reflect cumulative days
    X$time <- X$time + (i - 1) * 24
    
    # Store current day's results for the current dose
    dose_results <- rbind(dose_results, X)
  }
  
  # Add dose information to the results for faceting
  dose_results$Dose <- as.factor(Dose)
  
  # Combine current dose results with all results
  all_results <- rbind(all_results, dose_results)
}

# Convert days to hours for filtering
filtered_results <- subset(all_results, time <= 7 * 24)  

#
# NOTE: Below there is a line commented out that represents my attempt at overlaying 
# the observed data on top of the graph... it didn't work and needs a column rename 
# before running...
#

# plot linear y axis
ggplot(data = filtered_results, aes(x = time / 24, y = C1, color = "Concentration in C1")) + 
  geom_line() +
  #geom_point(data = mad, aes(x = time, y = C1, color = "Observed Concentration in C1")) +
  labs(title = "Concentration Over 7 Days for Multiple Oral Doses (Linear Scale)",
       x = "Time (days)", y = "Concentration (ng/L)", color = "Legend") +
  scale_x_continuous(breaks = seq(0, 7, by = 1),  # Set x-axis breaks for each day
                     labels = seq(0, 7, by = 1)) +   # Label breaks as day numbers
  scale_y_continuous(labels = label_comma()) +  # Use linear scale with normal notation
  facet_wrap(~ Dose, ncol = 1, labeller = labeller(Dose = label_both)) + # Add dose labels
  theme_minimal() +
  scale_color_manual(values = c("Simulated Concentration in C1" = "blue", 
                                "Observed Concentration in C1" = "red"))  # Customize colors
#
# NOTE: Below there is a line commented out that represents my attempt at overlaying 
# the observed data on top of the graph... it didn't work and needs a column rename 
# before running...
#

# plot log10 y axis
ggplot() +
  geom_line(data = filtered_results, aes(x = time / 24, y = C1, color = "Simulated Concentration in C1")) +
  # Overlay the observed data from mad, with points
  #geom_point(data = mad, aes(x = time, y = C1, color = "Observed Concentration in C1")) +
  labs(title = "Concentration Over 7 Days for Multiple Oral Doses (Log Scale)",
       x = "Time (days)", y = "Concentration (ng/L)", color = "Legend") +
  scale_x_continuous(breaks = seq(0, 7, by = 1),  # Set x-axis breaks for each day
                     labels = seq(0, 7, by = 1)) +   # Label breaks as day numbers
  scale_y_log10(labels = scales::label_comma()) +  # Use log scale with normal notation
  facet_wrap(~ Dose, ncol = 1, labeller = labeller(Dose = label_both)) + # Add dose labels
  theme_minimal() +
  scale_color_manual(values = c("Simulated Concentration in C1" = "blue", 
                                "Observed Concentration in C1" = "red"))  # Customize colors


#Part 6:
##========#========#========#========#========#========#========#========
#========Simulating the Pharmacodynamic response.#==========#============
##========#========#========#========#========#========#========#========

# Step 1: Observe the data
# Load data
mad_PD <- read.csv("cureaprobex_MAD_PD.csv")
head(mad_PD)  
Dose <- 0.1
mad_PD <- subset(mad_PD, dose == Dose )  

# Define the model with Emax effect on Wrkin1
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {
    
    dA_depot_dt <- -Ka * a_depot
    dC1_dt <- (Ka * a_depot - k12 * C1 + k21 * C2 - Ke * C1) / V1
    dC2_dt <- k12 * C1 - k21 * C2
    
    # show relationship between  emax and C1
    Effect <- Emax * C1 / (EC50 + C1)  # This reflects the effect on production, no additive 1 here
    
    dWrkin1_dt <- -k_inhibition * Wrkin1 * Effect
    
    # Return list of derivatives
    list(c(dA_depot_dt, dC1_dt, dC2_dt, dWrkin1_dt))
  })
}

# Define parameters for the model using Emax
theta <- list(
  Ka = fit$par[1],          
  V1 = fit$par[2],          
  k12 = fit$par[3],        
  k21 = fit$par[4],      
  Ke = fit$par[5],         
  Emax = 0.3,               # maximum effect 
  EC50 = 0.015,               # concentration at which half-maximal effect is achieved
  k_inhibition = 0.1        # rate of inhibition of Wrkin1 
)

#initial conditions
inits <- c(a_depot = Dose, 
           C1 = 0, 
           C2 = 0,
           Wrkin1 = 3) #I set 3 baserd on the observed data

#time settings
days <- 7         
times <- seq(0, 24 * days, by = 1) 

# run sim
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
  geom_line(data = mad_PD, aes(x = time_hr, y = wrkin1_pg_ml), color = "purple", size = 1, alpha = 0.7) +
  geom_point(data = mad_PD, aes(x = time_hr, y = wrkin1_pg_ml), color = "blue", size = 2, alpha = 0.7) +
  
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


# Part 7:
##========#========#========#========#========#========#========#========
#========#==========#========Dose recommendation.#==========#============
##========#========#========#========#========#========#========#========

# Step 1: Calculate the average of observed Wrkin1 concentration ignoring first data point and data before 4 hours
# filter  the observed data based on the dose and exclude first data points before a specific hour (outliers)
observed_wrkin1_values <- mad_PD$wrkin1_pg_ml[mad_PD$dose == Dose & mad_PD$time_hr >= 4][-1]

# observed data points average for 50%
average_observed_wrkin1 <- mean(observed_wrkin1_values, na.rm = TRUE)
dose_label <- paste("Simulated Response at", Dose, "mg")

# Create the plot overlaying observed and simulated data
ggplot() +
  # Observed data for the first 24 hours
  geom_line(data = mad_PD, aes(x = time_hr, y = wrkin1_pg_ml, color = "Observed Data"), size = 1, alpha = 0.7) +
  geom_point(data = mad_PD, aes(x = time_hr, y = wrkin1_pg_ml, color = "Observed Data"), size = 2, alpha = 0.7) +
  
  # Simulated Wrkin1 response for 24 hours
  geom_line(data = results_df, aes(x = time, y = Wrkin1, color = paste("Simulated Response at", Dose, "mg")), size = 1, linetype = "dashed") +
  
  #addline for 50% observed Wrkin1
  geom_hline(aes(yintercept = average_observed_wrkin1, color = "50% baseline Wrkin1"), linetype = "solid", size = 1.5) +
  
  # Customize the color scale for the legend
  scale_color_manual(name = "Legend",
                     values = c("Observed Data" = "black", 
                                "Simulated Response at 0.1 mg" = "blue",  # I need to change the label based on the dose
                                "50% baseline Wrkin1" = "red")) +
  
  # Theme and customization
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 24 * days, by = 6)) +  
  theme(legend.position = "top")  