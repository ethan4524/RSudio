library(tidyverse)
library(deSolve)
library(ggpubr)

rm(list = ls())

##Notes:
# Suproblatide: Precursor molecule that produces PRoblatide from the action of Extrase enzyme
# Problatide: A growth factor that stimulates cell division when bound to its receptor GF-R.
# Extrase: The enzyme responsible for converting Suproblatide to Problatide.
# GF-R:  The receptor that, when bound to Problatide, triggers cell division

#Specify deSolve Model
model <- function(time, inits, theta) {
  #Inputs:
  #time: vector of times to be included in output
  #y --> list of initial conditions
  #parms --> list of model parameters
  
  with(as.list(c(inits,theta)), {  #"with" get the named contents of y and theta and lets them be used inside the {}  
    
    #assume S >> ES,  S << Km,   S ~ St
    #compare a fresh file to this to find what we did in class
    ES = Et*St/(Km *(1 + Ie/Kie) + St)
    E = Et - ES
    S = St - ES
    fb = P / (Kd + P)
    
    dEt = (Prode * (fbss/fb)^n )- Kdege*E
    dSt = Prods - Kdegs*S - Kcat*ES
    dP = Kcat*ES - Kdegp * P
    
    #return list of derivatives
    list( c(dEt,dSt,dP))
  })
}

###################################################3
#Known constants
Km = 1.25 *1000 #pM
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
          n=2, #feedack
          Kd = Kd) #pM Product receptor binding affinity


# Simulation time (3 days)
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