library(tidyverse)
library(deSolve)
library(ggpubr)

rm(list = ls())

##Notes:
# Suproblatide: Precursor molecule that produces PRoblatide from the action of Extrase enzyme
# Problatide: A growth factor that stimulates cell division when bound to its receptor GF-R.
# Extrase: The enzyme responsible for converting Suproblatide to Problatide.
# GF-R:  The receptor that, when bound to Problatide, triggers cell division

# specify deSolve Model
model <- function(time, inits, theta) {
  with(as.list(c(inits, theta)), {  # Get named contents of inits and theta
    
    # Calculate enzyme-substrate complex and free enzyme
    ES = Et * St / (Km * (1 + Ie / Kie) + St)
    E = Et - ES
    S = St - ES
    
    # Fraction bound with inhibitor
    fb = (P) / (Kd + P)

    # Differential equations
    dEt = (Prode * (fbss / fb)^n) - Kdege * E
    dSt = Prods - Kdegs * S - Kcat * ES
    dP = Kcat * ES - Kdegp * P
    
    # Return list of derivatives
    list(c(dEt, dSt, dP))
  })
}

###################################################3
#Known constants
Km = 1.25 * 1000#nM 
Kcat = 4 #/sec
Kdege = 1 #/sec
Kdegs = 0.1 #/sec (the substrate can only be eliminated by converting into product)
Kd = 10 #pmol/L

#Concentrations measured at steady state
Pss = 20 #pmol/L
Sss = 10 #pmol/L
Etot_ss = 0.5 #pmol/L

#Calculated concentrations at SS

#Enzyme substrate complex at steady-state
ES_ss = Etot_ss * Sss/ (Km + Sss)

#Free enzyme at steady state
E_ss = Etot_ss - ES_ss

#Calculated production and elimination rates 
Prode = Kdege*E_ss    #Must equal free-enzyme elimination rate
Prods = Kdegs*Sss + Kcat*ES_ss  #Must equal free enzyme elimination rate + productformation rate
Kdegp = Kcat*ES_ss / Pss  #Calculate from steady-state condition for Product

fbss = Pss / (Pss + Kd)


#####################################3
inits = c(Et = Etot_ss,
          St = Sss  + ES_ss,
          P = Pss)

theta = c(Km = Km,
          Kcat = Kcat,
          Kdege = Kdege,
          Kdegs = Kdegs,
          Prode = Prode,
          Prods = Prods,
          Kdegp = Kdegp,
          Kie = 100, #pM
          Ie = 9500,
          Pss = Pss,
          n=3, #feedack
          Kd = Kd #pM Product receptor binding affinity
          )


# Simulation time
tlast = 3 * 24 * 60  # 3 days in minutes
times <- seq(0, tlast, 60)

#Simulate
X <- data.frame(ode(inits, times, model, theta))


#Calcualte fraction bound
X$fb = X$P / (theta["Kd"]  +X$P)



#Post-hoc calculations

# Plot Et (Extrase concentration)
plot_Et <- ggplot(X, aes(x = time, y = Et)) +
  geom_line(color = "blue",size = 1.2) +
  labs(x = "Time (min)", y = "Extrase Concentration (pmol/L)", title = "Extrase Concentration")

# Plot St (Suproblatide concentration)
plot_St <- ggplot(X, aes(x = time, y = St)) +
  geom_line(color = "orange",size = 1.2) +
  labs(x = "Time (min)", y = "Suproblatide Concentration (pmol/L)", title = "Suproblatide Concentration")

# Plot P (Problatide concentration)
plot_P <- ggplot(X, aes(x = time, y = P)) +
  geom_line(color = "red",size = 1.2) +
  labs(x = "Time (min)", y = "Problatide Concentration (pmol/L)", title = "Problatide Concentration")

# Plot Fraction Bound
plot_fb <- ggplot(X, aes(x = time, y = fb)) +
  geom_line(color = "purple") +
  labs(x = "Time (min)", y = "Fraction Bound", title = "Fraction Bound Over Time")

# Plot Fraction Bound vs Extrase Concentration
plot_fb_vs_Et <- ggplot(X, aes(x = fb, y = Et)) +
  geom_point(color = "blue") +
  labs(x = "Fraction Bound", y = "Extrase Concentration (pmol/L)", title = "Fraction Bound vs Extrase Concentration")

# Combine the first four plots into one grid and the fifth plot as a full plot below
ggarrange(plot_Et, plot_St, plot_P, plot_fb, 
          ncol = 2, nrow = 2, 
          common.legend = TRUE, legend = "bottom") %>%
  ggarrange(plot_fb_vs_Et, ncol = 1, nrow = 1, heights = c(2, 1))



#==============================Case 2, Part C==========================================#
# Create a vector of Fixaprob concentrations (pM)
I_values <- seq(0, 72, by = 1)  # Adjust as needed for the range you want

# Create an empty list to store the results
fraction_bound_results <- list()

# Loop through each Fixaprob concentration and run the model
for (I_conc in I_values) {
  
  # Update the inhibitor concentration (Kie) in the model parameters
  theta["Kie"] <- I_conc
  
  # Run the ODE simulation for this Fixaprob concentration
  X <- data.frame(ode(inits, times, model, theta))
  
  # Calculate fraction bound for each time point
  X$fb = X$P / (theta["Kd"] + X$P)
  
  # Store the time and fraction bound for this concentration
  fraction_bound_results[[as.character(I_conc)]] <- data.frame(Time = X$time, FractionBound = X$fb)
}

# Now convert the list to a data frame for easier plotting
fb_data <- do.call(rbind, lapply(names(fraction_bound_results), function(I_conc) {
  df <- fraction_bound_results[[I_conc]]
  df$FixaprobConcentration <- rep(I_conc, nrow(df))  # Add Fixaprob concentration column
  return(df)
}))

# Calculate the average fraction bound over time for each concentration
average_fb_data <- fb_data %>%
  group_by(FixaprobConcentration) %>%
  summarise(AverageFractionBound = mean(FractionBound, na.rm = TRUE))

# Plot fraction bound vs Fixaprob concentration (averaged over time)
ggplot(average_fb_data, aes(x = as.numeric(FixaprobConcentration), y = AverageFractionBound)) +
  geom_line(color = "blue", size = 1.2) +
  geom_vline(xintercept = 15, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0.222, color = "red", linetype = "dashed", size = 1) +
  labs(
    x = "Fixaprob Concentration (pM)",
    y = "Average Fraction Bound",
    title = "Average Fraction Bound vs Fixaprob Concentration"
  ) +
  annotate("text", x = 60, y = 0.25, label = "Target Fb = 0.222", color = "red", vjust = -0.5) +
  theme_minimal()
