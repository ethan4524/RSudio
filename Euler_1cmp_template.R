rm(list = ls())

library(tidyverse)
library(ggpubr)
#setwd("C:/Users/ethan/OneDrive/Documents/RStudio/RStudio")


#Define parameters
I =  5 #mg/hr
V = 5 #L
Ke = 10 #/hr
CL =  Ke*V  #L/hr
  
#Define Initial condition
C0 = 0

#Time step 
dt = 0.5  #hr

#Duration of simulation
tlast = 1  #hr

#Number of iterations
iterations = tlast / dt

#Empty variable to hold our results
xall = NULL

#Start with our state variable equal
#to the initial condition
C = C0
time = 0

for (i in 1:iterations) {
  #Save the current value of Concentration and time
  xall = rbind(xall, data.frame(Conc = C, Time = time)  )
  
  
  #Calculate the slope at the current 
  #timestep
  dcdt = (I - Ke*C)/V    #she said I/V - Ke*C
  
  #Calculate C at teh next timestep
  #based on the current value of C
  C = C + dcdt*dt
  
  #increase hte value of time to be used in the next timestep
  time = time + dt 
  print(C)
}

#Plot Results
ggplot(xall) + geom_path(aes(x=Time, y = Conc))



#Calculate analytical solution
t_analytical = seq(0,1,0.05)  #times at which to calculate the analytical solution.
C_analytical = (I/(V * Ke))*(1-exp(-Ke*t_analytical))
  
  
  
#Store analytical results in a dataframe. To make life easier later, we are using the same 
#column names that we used in the dataframe for our numerical results
x_analytical = data.frame(Time = t_analytical, Conc = C_analytical)


#Do a little work to get the numerical and analytical solution in the same dataframe
#Create a new column assigning a name to each simulation type:
xall$type = "Numerical"
x_analytical$type = "Analytical"

#Since both dataframes have the same column names, we can bind them together:
dat = rbind(xall, x_analytical)

#Now, plotting both is straightforward:

ggplot(dat) + geom_path(aes(x= Time, y = Conc, color = type))
  
dat1 = data.frame(time = time, data = xall, type = "Euler")
dat2 = data.frame(time = time, data = C, type = "Analytical")
dat = rbind(dat1, dat2)

ggplot(dat) + geom_path(aes ( x=time, y = data, color = type  ))

