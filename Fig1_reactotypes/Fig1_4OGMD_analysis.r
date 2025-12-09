library(dplyr)
library(magrittr)
library(readr)

paired_distances_df_m <- read_csv("~/.../paired_distances_df.m_Dgroup_noFMT.csv")
quantile(paired_distances_df_m$Paired_Distance, 0.75)
# 75% 
# 0.5716792 

paired_distances_df_m$distance_group <- ifelse(paired_distances_df_m$Paired_Distance >= 0.5716792, 
                                               "high_distance", 
                                               "low_distance")
library(dplyr)

# filter out distance_group_withFMT col and distance_group col, diff ones
paired_distances_df_m %>%
  filter(distance_group_withFMT != distance_group)

# Gut_Sample_ID Paired_Distance Oral_Sample_ID   Age Gender   BMI Group  MELD Severity_group reactotype distance_group_withFMT severity_group2
# <chr>                   <dbl> <chr>          <dbl> <chr>  <dbl> <chr> <dbl> <chr>               <dbl> <chr>                  <chr>          
#   1 RS_07_d01_ST            0.581 RS_07_d01_SA      70 M       36.6 LC       15 Low Severity            2 low_distance           >=15           
# 2 TW_01_ST                0.577 TW_01_SA          71 M       26.5 H        NA Healthy                 1 low_distance           NA  

write.csv(paired_distances_df_m,"~/.../paired_distances_df.m_Dgroup_noFMT2.csv",row.names=F)

################################################
############## OGMD vs MELD ####################
################################################

############## Linear model ####################

library(readr)
paired_distances_df.m <- read_csv("paired_distances_df.m_Dgroup_noFMT2.csv")

### linear MELD ~ oral-gut_distance
# Create the scatter plot
ggplot(paired_distances_df.m, aes(x = MELD, y = Paired_Distance)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Column 1", y = "Column 2", title = "Linear Correlation Plot")+
  theme_minimal()

# Using Pearson correlation:
cor.test(paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),]$MELD, 
         paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),]$Paired_Distance, method = "pearson")
# Using Spearman correlation:
cor.test(paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),]$MELD, 
         paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),]$Paired_Distance, method = "spearman")

### linear model
lm(MELD ~ Paired_Distance, data = paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),])
model <- lm(Paired_Distance ~ MELD, data = paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),])
summary(model)
model_1 <- glm(Paired_Distance ~ MELD, data = paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),])
summary(model_1)

# visualisation
# Create a scatter plot with linear regression line
ggplot(paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),], aes(x = MELD , y = Paired_Distance)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = T) +  # Add linear regression line
  labs(x = "MELD", y = "Oral-gut Reactobiome Distance") +  # Add axis labels
  ggtitle("Linear Regression between MELD and OGD") +  # Add title
  theme_bw()  # Use a minimal theme

# First, make sure you have the ggpmisc package installed
# install.packages("ggpmisc")

library(ggplot2)
library(ggpmisc)

# Creating the plot with p-value and R² displayed


p1 <- ggplot(paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),], aes(x = MELD , y = Paired_Distance)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = TRUE) +  # Add linear regression line
  labs(x = "MELD", y = "Oral-gut Reactobiome Distance") +  # Add axis labels
  ggtitle("Linear Regression between MELD and OGD") +  # Add title
  theme_bw() +  # Use a minimal theme
  stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
               label.x.npc = 0.8, label.y.npc = 0.2,  # Adjust these values for the p-value label position
               formula = y ~ x, parse = TRUE) 

p1 # save Linear Regression between MELD and OGD 5.5*5 L

monet_palette <- c("Average_High" = "#809FFF",  # LightSkyBlue
                   "Average_Low" = "#00008B")  # LightPink

# color by distance_group

ggplot(paired_distances_df.m[-which(is.na(paired_distances_df.m$MELD)),], 
       aes(x = MELD, y = Paired_Distance)) +  
  geom_point(aes(color = distance_group), size = 4) +  # Color points by distance_group and set size 点的大小
  geom_smooth(method = "lm", se = TRUE, color = "black") +  # Add single linear regression line
  labs(x = "MELD", y = "Oral-gut Reactobiome Distance") +  # Add axis labels
  ggtitle("Linear Regression between MELD and OGD") +  # Add title
  theme_bw() +  # Use a minimal theme
  scale_color_manual(values = c("high_distance" = "#809FFF", "low_distance" = "#00008B")) +  # Specify colors
  stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
               label.x = 30, label.y = 2,  # Adjust these values for the p-value label position based on data coordinates
               formula = y ~ x, parse = TRUE)


