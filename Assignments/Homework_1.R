rm(list = ls())

library(tidyverse)
library(ggpubr)
setwd("C:/Users/Ethan/OneDrive/Documents/BIOE8510/RSudio/Assignments")

#-------------PART 1----------------#
####Problem 1
#A: [Enter description here]
help("sort")
#B: [Enter description here]
help("rnorm")
#C:
myRandomNums=rnorm(100,mean=50,sd=25)
#D:
head(myRandomNums)
tail(myRandomNums)
#E:
mySortedSums=sort(x=myRandomNums,decreasing=TRUE)
#F:
length(mySortedSums)
#G:
plot(mySortedSums)
#H:
hist(mySortedSums)

####Problem 2
a=110^(1/3)
b=tan(45 * (pi/180))
c=pi^2
d=log(20.1)
e=3.3^4
A=list(a,b,c,d,e)

####Problem 3
# A
v=seq(0,100,by=10)
time = c(v,v,v,v)
# B
expnum=c()
for (i in 1:4) {
  for (j in 0:10) {
    expnum=c(expnum,i)
  }
}
# C
baseline=runif(length(time))
# D
mydf=data.frame(time,expnum,baseline)
# E
print(colnames(mydf))
# F
dim(mydf)
# G
response = rnorm(n = length(mydf$time), mean = 5, sd = 0.5)
mydf['response'] = response
# H
mydf['difference'] = mydf$response - mydf$baseline
# I
exp2 = subset(mydf, expnum==2)
# J
time0 = subset(mydf, time==0)
# K
exp4=subset(mydf,expnum==4)
temp=data.frame(subset(exp4,time==50))
print(temp$difference)
# L
responders=subset(subset(mydf,response>5), time==100)
# M
testsample=subset(mydf, (response < 5) & (baseline > 0.5))
# N
mean_difference=mean(mydf$difference)
# O
sum_baseline=sum(mydf$baseline)
# P
new_row = data.frame(
  time=0,
  expnum=5,
  baseline=0,
  response=5,
  difference=5
)
mydf = rbind(mydf, new_row)
# Q
exp1_70=(subset(mydf, (expnum==1 & time==70)))$response
exp3_70=(subset(mydf, (expnum==3 & time==70)))$response
diff=(exp1_70-exp3_70)

####Problem 4
#A
Loblolly
#B
head(Loblolly)
#C
unique_seed_sources=length(unique(Loblolly$Seed))
#D (note: there was a strange issue with the factor being preserved. Therefore, I converted the value into a numeric)
age_10=subset(Loblolly, age==10)
row = age_10[which.min(age_10$height),]
seed_10=as.numeric(as.character(row$Seed))
#E
S301 = subset(Loblolly, Seed==301 )
#F
S323 = subset(Loblolly, Seed==323 )
#G
LoblollySmall = rbind(S301,S323)

#-------------PART 2----------------#
# Notes: Use ggplot for creating plots
# For saving images, click export over the plot in the plot window
# ggarrage(plot1, plot2) will combine two plots (from the gg package)

#A  
# Is it easy to see trends in the data that might help you distinguish
# carcinoma from other types?

bc = read.csv("breastcancer.csv")

#B
s1=ggplot() + geom_point(data=bc,mapping=aes(x=PA500, y=I0))
#C
s1=ggplot() + geom_point(data=bc,mapping=aes(x=PA500, y=I0, color=Class))
#D
s1=ggplot() + geom_point(data=bc,mapping=aes(x=PA500, y=I0, color=Class, size=PA500))
# Entries with higher I0 tend to have a smaller PA500 value.
# The class adi and the class con have a much higher I0 value.
# Class car had the greatest values for PA500 and a similar I0 value compared to classes excluding adi and con

#E
# Repeating step D:
#x=HFS, y = I0
s2=ggplot() + geom_point(data=bc,mapping=aes(x=HFS, y=I0, color=Class))
#x=DA, y = I0 ; No interesting trends 
s3=ggplot() + geom_point(data=bc,mapping=aes(x=DA, y=I0, color=Class))
# when comparing DA to I0, the resulting graph separated each class into their
# own regions on the scatterplot based on I0. 
# class adi had the greatest values for I0 and the greatest spread for DA.
# mas, gla, and car all have similar values for da and I0.
s4=ggplot() + geom_point(data=bc,mapping=aes(x=I0, y=Max.IP, color=Class))
# Trend: I0 and Max.IP are proportional 

# Feature Classification
#A
s5=ggplot() + geom_point(data=bc,mapping=aes(x=PA500, y=I0, color=Class))
# Using the plot s5, its safe to assume that a
# class with a PA500 value >= 0.2 is carninoma.
# Classes below 0.2 PA500 value are non-carcinoma.

#B 
