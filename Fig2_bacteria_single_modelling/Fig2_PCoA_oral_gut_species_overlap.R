#### remove all variables
# rm(list = ls())
## load package ###
library(readxl)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(dplyr)
library(vegan)
library(tibble)

# PCoA Oral Gut Reactobiome
library(readr)
reactobiome_oral <- read_csv("~/.../reaction_relab_frommsp_oral_lc_and_h_noFMT.csv")
reactobiome_gut <- read_csv("~/.../reaction_relab_frommsp_gut_lc_and_h_noFMT.csv")

library(readxl)
Paired_LC_metadata <- read_excel("~/.../Paired_LC_metadata.xlsx")

oral_ID = Paired_LC_metadata$Oral_Sample_ID
gut_ID = Paired_LC_metadata$Gut_Sample_ID

reactobiome_oral.paired = reactobiome_oral[,c("rn",oral_ID)] # [1] 2772   54
reactobiome_gut.paired = reactobiome_gut[,c("rn",gut_ID)] # [1] 3608   54

reactobiome_OG = merge(reactobiome_oral.paired,reactobiome_gut.paired,by="rn",all=T) # [1] 3612  107
reactobiome_OG[is.na(reactobiome_OG)]=0

reactobiome_OG.n = data.matrix(reactobiome_OG[,-1])
rownames(reactobiome_OG.n) = reactobiome_OG$rn

reactobiome_OG.n.t = t(reactobiome_OG.n) %>% as.data.frame() # [1]  106 3612

metadata_OG = rbind(Paired_LC_metadata,Paired_LC_metadata) # [1] 106  19
metadata_OG$Sample_ID = rownames(reactobiome_OG.n.t)
metadata_OG$Site <- c(rep("Saliva", 53), rep("Feces", 53))

reactobiome_OG_df = cbind(metadata_OG,reactobiome_OG.n.t)

###### ---------------------------------------------------- #######
# df metadata
head(reactobiome_OG_df)

# data matrix
reactobiome_OG.n.t # rows are samples, and columns are features

###### ---------------------------------------------------- #######
######-------- MAKE PCoA plots ---------------------------- #######
distance <- vegdist(reactobiome_OG.n.t,method = "bray") # calculate the distance between samples (use method "bray")
pcoa <- cmdscale(distance,k=2,eig = TRUE) # returns the best-fitting k-dimensional representation, where k may be less than the argument k.(tips: return ?cmdscale to see more in Help Window)
plot_data <- data.frame({pcoa$points}) # Convert the points from the pcoa analysis in the previous step into a data frame and assign it to plot_data, which will be used when drawing the plot
eig <- pcoa$eig #Assign a value to the eig in the pcoa processing return, which is needed to calculate the confidence interval below
data <- data.frame(reactobiome_OG_df$Site,reactobiome_OG_df$Cohort,reactobiome_OG_df$distance_group,reactobiome_OG_df$Severity_group,plot_data) # Assign the two grouped columns and the pcoa-processed point data to data for graphing purposes
names(data) <- c("Site","Cohort","Distance_group","Severity_group","PCoA1",'PCoA2')

# order the levels
data$Site <- factor(data$Site, levels = c("Saliva", "Feces"))
levels(data$Site) # [1] "Saliva" "Feces"

data$Distance_group <- factor(data$Distance_group, levels = c("low_distance", "high_distance"))
levels(data$Distance_group) # [1] "low_distance"  "high_distance"

plot_PCoA <- 
  ggplot(data = data, aes(x = PCoA1, y = PCoA2, shape = Site, color = Distance_group)) +   
  # Base layer (scatter plot)
  geom_point(alpha = 1, size = 2) +    
  
  # Ellipses with dashed lines
  stat_ellipse(aes(color = Site), type = "norm", geom = "path", linetype = "dashed", size = 1) +    
  
  # Custom labels for axes, including PCoA variance explained
  labs(
    x = paste("PCoA1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
    y = paste("PCoA2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = "")
    # title = "HC"
  ) +    
  
  # Specify color for Distance_group (this stays the same)
  scale_color_manual(values = c("high_distance" = "#00008B", "low_distance" = "#809FFF")) +    
  
  # Remove the titles of legends
  guides(
    color = guide_legend(title = NULL),  # Remove title for color legend (Distance_group)
    shape = guide_legend(title = NULL)   # Remove title for shape legend (Site)
  ) +    
  
  # Custom theme for axis text, legend, and background
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(colour = 'black', size = 15),
    axis.title.y = element_text(colour = 'black', size = 15),
    legend.text = element_text(size = 12)
  )

plot_PCoA # save PCoA_OG_LC 4.5*3 P
















