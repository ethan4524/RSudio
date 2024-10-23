#Load packages
library(lattice)
library(tidyverse)
library(deSolve)



#Clear out the workspace
rm(list=ls())


#Specify deSolve Model
model <- function(time, inits, theta) {
  #Inputs:
  #time: vector of times to be included in output
  #y --> list of initial conditions
  #parms --> list of model parameters
   
  with(as.list(c(inits,theta)), {  #"with" get the named contents of y and theta and lets them be used inside the {}  
    
    CL = Ke*V

        
    #Model equation
    dC <- (I - CL*C)/V
    
    
    #return list of derivatives
    list( c(dC))
  })
}

#Set initial conditions
#IMPORTANT: the variable in your initial condition list must have exactly the same name as the variable in your model
inits <- c(C =0 )   

#Make a named list of model parameters
theta <- c(I = 5,   #mg/hr
           V = 5,   #L
           Ke = 10)  #/hr  #list of parameters

#Simulation time
tlast = 1 ; #hr

times <- seq(0, tlast, 0.01)

#Simulate
X <- ode(inits, times, model, theta)

#convert results to a dataframe to make plotting easier
X = data.frame(X)


#Plot 
dev.new()
ggplot(X) + geom_path(mapping = aes(x=time, y = C)) + xlab("time") + ylab("C")


###############################################################################
###############################################################################
#Plot numerical and analytical solution

xall = X


I = 5 #mg/hr
V = 5 #L
Ke = 10 #/hr
CL = Ke*V #L/hr



#Calculate analytical solution
t_analytical = seq(0,1,0.01)  #times at which to calculate the analytical solution.
C_analytical = (I/CL)*(1-exp(-(CL/V)*t_analytical))

#Store analytical results in a dataframe. To make life easier later, we are using the same 
#column names that we used in the dataframe for our numerical results
x_analytical = data.frame(time = t_analytical, C = C_analytical)

#Do a little work to get the numerical and analytical solution in the same dataframe
#Create a new column assigning a name to each simulation type:
xall$type = "Numerical"
x_analytical$type = "Analytical"

#Since both dataframes have the same column names, we can bind them together:
dat = rbind(xall, x_analytical)

#Now, plotting both is straightforward:

ggplot(dat) + geom_path(aes(x= time, y = C, color = type)) + 
  geom_point(aes(x= time, y = C, color = type))





###############################################################################
###############################################################################

# Infuse for 30 minutes and then turn off infusion


###### Infuse for 30 minutes

inits <- c(C =0 )   

#Make a named list of model parameters
theta <- c(I = 5,   #mg/hr
           V = 5,   #L
           Ke = 10)  #/hr  #list of parameters

#Simulation time - 30 minutes 
tlast = 0.5 ; #hr

times <- seq(0, tlast, 0.01)

#Simulate
X <- ode(inits, times, model, theta)

#convert results to a dataframe to make plotting easier
X = data.frame(X)


#####################Turn off infusion and simulate for 1.5 hrs

theta["I"] = 0

inits <- c(C =   tail(X$C, n=1)  ) 

times <- seq(0, 1.5, 0.01)

#Simulate
X2 <- ode(inits, times, model, theta)

X2 = data.frame(X2)

#Shift X2 to the right by half an hour
X2$time = X2$time + 0.5

ggplot(X2) + geom_path(mapping = aes(x=time, y = C)) + xlab("time") + ylab("C")


#Combine first 30 minutes and last 90 minutes into one data frame

xall = rbind(X, X2)
ggplot(xall) + geom_path(mapping = aes(x=time, y = C)) + xlab("time") + ylab("C")



