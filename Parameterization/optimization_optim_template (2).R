library(tidyverse)

setwd("C:/Users/hallowkm/Documents/Teaching/BIOE3720 S2021/Module 7")

#Step 1:
#Load data set
#t --> time, independent variable
#R --> response, dependent variable
A = read.table("practiceData.csv", header= T, sep = ",")

#Step 2: 
#Explore Data - what do the relationships between variables in your dataset look like? 


#Step 3: Specify a model. Usually start with the simplest model that will:
#1) adequately describe the relationship between variables, and
#2) allow you to pursue the purpose of your modeling activity (describe? Explain? Predict?)

#Model 1: linear  (we already know this model is too simple for our data)
b = 0#intercept initial guess
m = 15 #slope initial guess

#Specify model for R, and calculate values using initial guesses
#Calculate R at each time t for which R was measured
A$Rmod = 
  
#Plot model and data
ggplot(A) + geom_point(aes(x=t, y=R)) + 
		geom_path(aes(x=t, y = Rmod))


#Calculate error
MSE = 

print(MSE)
####Try Rerunning for different initial guesses above. How does it change the fit? How does the error change?


  
################### Optimization using Optim ########################

#Guess starting values for parameters:
theta = c(m=15, b = 0)


#Define an objective function. This function should take in a list of parameters as input and return the
#objective function value

objfxn <- function(theta) {
  #theta is a list of model parameters
  
  #Calculate model outputs for input parameters theta
  A$Rmod = 
  
  #Calculate objective function. In this case, we want to minimize the MSE, so our objective function value is the MSE
  obj = 
  
  #For debugging purposes, it can be useful to print the objective function and parameter values  
  print(objval)
  print(theta)
  
  #Function returns the objective function value
  return(obj)
}

#test objfxn

#Test your objective function by running with a test case. It should return a single value - the MSE.
objfxn(theta) #Running this line should return a single value for your objective function.



fit = optim(theta, objfxn, method = "BFGS", hessian = T)


#Results of the optimization are stored in the variable "fit"
fit



#Update theta with your optimized parameters
fit$par
theta$m = fit$par[1]
theta$b = fit$par[2]


#Check 1: did it converge?
fit$convergence   #If this is anything but zero, the optmization did not converge - it was not able to find the minimum value for the objective function

#Check 2: Visual Check. Plot observed and predicted data

#Calculate model output with final parameter values
A$Rmod_fit = #you fill in here ....
  


#Check 3: Residuals - The residuals are the differences between observed and predicted values. 

A$resid = #fill in here ....
  
  
  #histogram of residuals
  hist(A$resid)

#Residuals vs. independent variable (time)
# - add a horizontal reference line at 0



#Check 4: Observed vs. Predicted
# - add a diagonal reference line with a slope of 1


#Homework: Improve model structure, optimize again, and replot graphical checks. 
#Hint: Use logistic function instead of linear function:   R = Rmax/(1+exp(-k*(t-t50)))

#Rmax --> max R as t goes to infinity
#k --> growth rate constant
#t50 --> time of half maximum response


