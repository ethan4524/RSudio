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
