# Clear the workspace
rm(list = ls())

# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Set working directory
setwd("C:/Users/Ethan/Documents/BIOE8510/RSudio/Assignments")

# Read the data
bp = read.csv("BP.csv")

# a: Calculate the mean blood pressure for each ID
part_a <- bp %>% 
  group_by(ID) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE))

# b: Calculate the mean, standard deviation, and standard error for each CAT and DIET pair
part_b <- bp %>%
  group_by(CAT, DIET) %>%
  summarise(
    mean_bp = mean(value, na.rm = TRUE),
    sd_bp = sd(value, na.rm = TRUE),
    n = n(),
    sem_bp = sd_bp / sqrt(n)
  )
# Create the paired bar graph
plot4 <- ggplot(part_b, aes(x = CAT, y = mean_bp, fill = DIET)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Bars for each diet within each strain
  geom_errorbar(aes(ymin = mean_bp - sem_bp, ymax = mean_bp + sem_bp), 
                position = position_dodge(0.7), width = 0.25) +    # Error bars for SEM
  labs(title = "Mean Blood Pressure by Strain and Diet",
       x = "Strain",
       y = "Mean Blood Pressure (mmHg)") +
  theme_minimal()

# Display the plot
print(plot4)

# Optionally save the plot as a file
# ggsave(filename="plot4.png", plot=plot4, width=8, height=6, dpi=300)
