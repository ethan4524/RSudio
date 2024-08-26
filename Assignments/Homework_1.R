rm(list = ls())

library(tidyverse)
library(ggpubr)
setwd("C:/Users/Ethan/OneDrive/Documents/BIOE8510/Assignments")

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

