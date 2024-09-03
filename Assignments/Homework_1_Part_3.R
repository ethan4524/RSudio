rm(list = ls())

library(tidyverse)
library(ggpubr)
setwd("C:/Users/Ethan/OneDrive/Documents/BIOE8510/RSudio/Assignments")

bp = read.csv("BP.csv")
head(bp)
#Plot 1:
bp=subset(bp,time<0.8)

plot1 <- ggplot(bp, aes(x = time, y = value, color = DIET, group = ID)) +
  geom_line() +                            # Add lines for each ID
  stat_summary(fun = mean, geom = "line", aes(group = DIET,color=ID), color = "black", size = 1) +  # Optional: Add mean line per diet
  labs(title = "Blood Pressure Over Time by Diet and ID",
       x = "Time (seconds)",
       y = "Blood Pressure (mmHg)") +
  theme_minimal()
ggsave(filename="plot1.png",plot=plot1, width=8,height=6, dpi=300)
plot1

plot2=ggviolin(bp, x = "DIET", y = "value", fill = "CAT",
         ylab = "Blood Pressure (mmHg)", xlab = "Diet",
         title = "Blood Pressure Distribution by Diet and Strain",
         add = "boxplot", add.params = list(fill = "white"))
ggsave(filename="plot2.png",plot=plot2, width=8,height=6, dpi=300)
plot2
plot3=ggplot(bp, aes(x = time, y = value, group = ID, color = ID)) +
  geom_line() +
  facet_wrap(~ DIET) +
  labs(title = "Blood Pressure Over Time by Diet",
       x = "Time (seconds)",
       y = "Blood Pressure (mmHg)") +
  theme_minimal()
ggsave(filename="plot3.png",plot=plot3, width=8,height=6, dpi=300)
plot3


bp_summary <- aggregate(value ~ DIET + CAT, data = bp, FUN = mean)

# Bar plot of average blood pressure by Diet and Strain
plot4=ggbarplot(bp_summary, x = "DIET", y = "value", fill = "CAT",
          ylab = "Average Blood Pressure (mmHg)", xlab = "Diet",
          title = "Average Blood Pressure by Diet and Strain",
          position = position_dodge(0.8))
ggsave(filename="plot4.png",plot=plot4, width=8,height=6, dpi=300)
plot4

plot5=ggplot(bp, aes(x = time, y = value, group = ID, color = ID)) +
  geom_line() +
  facet_grid(CAT ~ DIET) +               # Facet grid by CAT (rows) and DIET (columns)
  labs(title = "Blood Pressure Trends Over Time by Diet and Strain",
       x = "Time (seconds)",
       y = "Blood Pressure (mmHg)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) 
ggsave(filename="plot5.png",plot=plot5, width=8,height=6, dpi=300)
plot5