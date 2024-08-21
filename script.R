library(tidyverse)
library(ggpubr)

setwd("C:/Users/Ethan/OneDrive/Documents/BIOE8510")

df = read.csv("breastcancer.csv")

head(df)

# Below is an example of tidyverse function %>% as well as the select and filter function
df_Class = df %>% select(Class) %>% filter(Class=="car")

print(df_Class)
