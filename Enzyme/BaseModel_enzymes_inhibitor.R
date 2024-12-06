#Load packages


library(tidyverse)
library(deSolve)
library(ggpubr)


#Specify deSolve Model
model <- function(time, inits, theta) {
  #Inputs:
  #time: vector of times to be included in output
  #y --> list of initial conditions
  #parms --> list of model parameters
  
  with(as.list(c(inits,theta)), {  #"with" get the named contents of y and theta and lets them be used inside the {}  
    
    #assume S >> ES,  S << Km,   S ~ St
    #compare a fresh file to this to find what we did in class
    ES = Et*St/(Km *(1 + Ie/Kie) + St)  #<--- This could be wrong but 99% sure not
    
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
Km = 0.75*1000 #pM
Kcat = 1 #/sec
Kdege = 1 #/sec
Kdegs = 1 #/sec (the substrate can only be eliminated by converting into product)
Kd = 20 #pM

#Concentrations measured at steady state
Pss = 20 #pmol/L
Sss = 10 #pmol/L
Etot_ss = 1 #pmol/L

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
          n=1, #feedack
          Kd = Kd) #pM Product receptor binding affinity



#Simulation time
tlast = 30*60*24 ; 
times <- seq(0, tlast, 60)

#Simulate
X <- data.frame(ode(inits, times, model, theta))


#Calcualte fraction bound
X$fb = X$P / (theta["Kd"]  +X$P)

#Post-hoc calculations

#Plot
ggplot(X) + geom_path(aes(x=time, y = Et))
ggplot(X) + geom_path(aes(x=time, y = St))
ggplot(X) + geom_path(aes(x=time, y = fb))

ggplot(X) + geom_path(aes(x=time, y = P))

tail(X$P)


#####################################
#Play around with changing parameters
#theta["Prods"] = Prods


#Adding an inhibitor
theta["I"] =1000

theta["n"] = 1


#Simulate
X <- data.frame(ode(inits, times, model, theta))


#Post-hoc calculations

#Plot
ggplot(X) + geom_path(aes(x=time, y = P)) + geom_hline(yintercept = 2, linetype = "dashed")

#######################################################################################





