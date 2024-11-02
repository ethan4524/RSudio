#Load packages
library(lattice)
library(ggplot2)
library(deSolve)

setwd("C:/Users/ethan/OneDrive/Documents/RStudio/RSudio/In Class/Pharmacokdynamics")

HR = read.csv("Theoph_PD_HR.csv")

ggplot(HR) + geom_point(aes(x=Time_hours, 
                            y=HR_bpm, 
                            color = factor(Dose_mg), 
                            group = Dose_mg))
#HR = HR0*(1+ Emax*Concentration) / (Km + Concentration)
#Km = EC50 = IC50 (this is is just the concentraion that gives you half of the max response)

########################################### Define 1 compartment PK Model with oral dosing #####################

#Specify deSolve Model
model <- function(time, inits, theta) {
  #Inputs:
  #time: vector of times to be included in output
  #y --> list of initial conditions
  #parms --> list of model parameters
  
  with(as.list(c(inits,theta)), {  #"with" get the named contents of y and theta and lets them be used inside the {}  
    
    Effect = 1+Emax * C1 / (EC50 + C1)
    
    #Model equation
    dA_depot = -Ka*A_depot
    dC1 <- Ka*A_depot/V1 - Ke*C1 
    
    #Add the dif eq for PD model HR into this model instead of a separate model:
    #Effect = 1+Emax*C/(C+Km)
    
    dHR <- Kin*Effect-Kout*HR
    
    #return list of derivatives
    list( c(dA_depot, dC1,dHR))
  })
}

#Fit model parameters
theta = c(Ka = 1.777, V1 =29.39, Ke = 0.054)

#Selecting Subject 1
Theoph_1 = subset(Theoph,Subject==1)

#Note that in Theophylline dataset, Dose column is given as mg/kg. Total dose is Dose*Wt
#Calculate a representative dose to use while exploring:
#Dose = Theoph_1$Dose[1]*Theoph_1$Wt[1]
Dose = 1000
#Inits
inits <- c(A_depot = Dose, #the amount in our depot compartment a time zero is the dose for oral dosing
           C1 =   0 ,
           HR = 70
           )    
#Simulation time
tlast = 24 ; 
times <- seq(0, tlast, 0.01)


theta["Emax"] = 0.5
theta["EC50"] = 4.8
theta["HR0"] = 70
theta["Kin"] =42
theta["Kout"] = theta["Kin"] / inits["HR"] #by steady state these are non independant parameters 

#Simulate
X <- ode(inits, times, model, theta)

#convert results to a dataframe to make plotting easier
X = data.frame(X)


#HR = HR0*(1+ (Emax*Concentration) / (Km + Concentration))
X$HR = theta["HR0"] * (1 + (theta["Emax"] * X$C1)/(theta["EC50"]  + X$C1 ))

ggplot() + geom_point(data = subset(HR,Dose_mg==1000), 
                      aes(x=Time_hours,y=HR_bpm,color = "red")) +
  geom_path(data = X, aes(x=time,y=HR))

# Part 2:

doses = unique(HR$Dose_mg)

Xall = NULL

for (i in 1:length(doses)) {
  
  #Simulate
  X <- ode(inits, times, model, theta)
  
  #convert results to a dataframe to make plotting easier
  X = data.frame(X)
  X$dose = doses[i]
  
  Xall = rbind(Xall,X)
  
}


#Plot simulation on top of Theophylline data:
#note

ggplot() + geom_point(data = subset(HR,Dose_mg==Dose), 
                      aes(x=Time_hours,y=HR_bpm,color = "red")) +
  geom_path(data = X, aes(x=time,y=HR)) + 
  ylim(c(65,105))

ggplot() + geom_point(data = HR, aes(x=Time_hours, y=HR_bpm,color=factor(Dose_mg))) +
  geom_path(data = Xall, aes(x=time,y=HR,color=factor(dose)))
  

                   


