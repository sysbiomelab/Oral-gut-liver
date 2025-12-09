######################################################
####### OGMD differential abundance reactions ########
######################################################


######### load relative reaction abundance #######
library(readr)
reaction_relab_frommsp_gut_lc_and_h <- read_csv("~/.../reaction_relab_frommsp_gut_lc_and_h_noFMT.csv")
reaction_relab_frommsp_oral_lc_and_h <- read_csv("~/.../reaction_relab_frommsp_oral_lc_and_h_noFMT.csv")

######### sig high vs low oral-gut-distance #######
#### make sample id index ####
pair_oral_high_id = paired_distances_df_m[which((paired_distances_df_m$distance_group == "high_distance")&
                                                  (paired_distances_df_m$Group == "LC")),]$Oral_Sample_ID
pair_oral_low_id = paired_distances_df_m[which((paired_distances_df_m$distance_group == "low_distance")&
                                                 (paired_distances_df_m$Group == "LC")),]$Oral_Sample_ID

pair_gut_high_id = paired_distances_df_m[which((paired_distances_df_m$distance_group == "high_distance")&
                                                 (paired_distances_df_m$Group == "LC")),]$Gut_Sample_ID
pair_gut_low_id = paired_distances_df_m[which((paired_distances_df_m$distance_group == "low_distance")&
                                                (paired_distances_df_m$Group == "LC")),]$Gut_Sample_ID

pair_oral_id = paired_distances_df_m[which(paired_distances_df_m$Group == "LC"),]$Oral_Sample_ID

pair_gut_id = paired_distances_df_m[which(paired_distances_df_m$Group == "LC"),]$Gut_Sample_ID

#### make data matrix #####
gut_rxn = reaction_relab_frommsp_gut_lc_and_h[,-1] %>% t() %>% data.matrix()
colnames(gut_rxn) = reaction_relab_frommsp_gut_lc_and_h$rn
gut_rxn[,]

oral_rxn = reaction_relab_frommsp_oral_lc_and_h[,-1] %>% t() %>% data.matrix()
colnames(oral_rxn) = reaction_relab_frommsp_oral_lc_and_h$rn

#### sig test ####
#### wilcoxon test for Reaction ids ####
### oral ###
pair_oral_high_id
pair_oral_low_id
pair_oral_id
oral_rxn
oral_paired_rxn = oral_rxn[pair_oral_id,] %>% data.matrix()

### wilcoxon
# pvalue
p_values <- vector("list", ncol(oral_paired_rxn))
for(i in seq_along(1: ncol(oral_paired_rxn))){
  p_values[i] =  wilcox.test(oral_paired_rxn[pair_oral_high_id,i], oral_paired_rxn[pair_oral_low_id,i],paired=F,exact=FALSE) $p.value}
p_values = data.frame(p_values = sapply(p_values, c))
p_values$rxn = colnames(oral_paired_rxn)
p_values$FDR = p.adjust(p_values$p_values)

p_values_oral = p_values

# Log2FC
high_mean <- numeric(length = ncol(oral_paired_rxn))
low_mean <- numeric(length = ncol(oral_paired_rxn))
Log2FC <- numeric(length = ncol(oral_paired_rxn))

for (i in seq_along(1:ncol(oral_paired_rxn))) {
  high_mean[i] = mean(oral_paired_rxn[pair_oral_high_id, i])
  low_mean[i] = mean(oral_paired_rxn[pair_oral_low_id, i])
  Log2FC[i] = log((low_mean[i] + 0.000001) / (high_mean[i] + 0.000001), 2)
}

FC1.1 = data.frame(high_mean = sapply(high_mean, c),
                   low_mean = sapply(low_mean, c),
                   Log2FC = sapply(Log2FC, c))

oral_rxn_stat = cbind(p_values_oral,FC1.1)
head(oral_rxn_stat)

#### wilcoxon test for Reaction ids ####
###  gut #####
pair_gut_high_id
pair_gut_low_id
pair_gut_id
gut_rxn
gut_paired_rxn = gut_rxn[pair_gut_id,] %>% data.matrix()

### wilcoxon
# pvalue
p_values <- vector("list", ncol(gut_paired_rxn))
for(i in seq_along(1: ncol(gut_paired_rxn))){
  p_values[i] =  wilcox.test(gut_paired_rxn[pair_gut_high_id,i], gut_paired_rxn[pair_gut_low_id,i],paired=F,exact=FALSE) $p.value}
p_values = data.frame(p_values = sapply(p_values, c))
p_values$rxn = colnames(gut_paired_rxn)
p_values$FDR = p.adjust(p_values$p_values)

p_values_gut = p_values
head(p_values_gut)

# Log2FC
high_mean <- numeric(length = ncol(gut_paired_rxn))
low_mean <- numeric(length = ncol(gut_paired_rxn))
Log2FC <- numeric(length = ncol(gut_paired_rxn))

for (i in seq_along(1:ncol(gut_paired_rxn))) {
  high_mean[i] = mean(gut_paired_rxn[pair_gut_high_id, i])
  low_mean[i] = mean(gut_paired_rxn[pair_gut_low_id, i])
  Log2FC[i] = log((low_mean[i] + 0.000001) / (high_mean[i] + 0.000001), 2)
}

FC1.2 = data.frame(high_mean = sapply(high_mean, c),
                   low_mean = sapply(low_mean, c),
                   Log2FC = sapply(Log2FC, c))

gut_rxn_stat = cbind(p_values_gut,FC1.2)
head(gut_rxn_stat)

#### save table ###
# setwd("~/Documents/github_reactobiome/reactobiome/rytoPJ_231013/reactobiome_231013/O-G_distance_stat_R/output_tables")
write.csv(gut_rxn_stat,
          "~/.../gut_Rid_stat_no_health.csv",
          row.names = FALSE)
write.csv(oral_rxn_stat,
          "~/.../oral_Rid_stat_no_health.csv",
          row.names = FALSE)


######################################################
####### plot differential reactions ##################
######################################################

library(readxl)
OGD_low_high_stat_oral <- read_excel("~/.../Oral_gut_Rid_Low_vs_High.xlsx", 
                                     sheet = "oral")
OGD_low_high_stat_gut <- read_excel("~/.../Oral_gut_Rid_Low_vs_High.xlsx", 
                                    sheet = "gut")

### use not cln one, and add a row cbind rn and Rid ####
OGD_low_high_stat_oral <- OGD_low_high_stat_oral %>%
  mutate(rxnInfo = paste(rxn, rn, sep = " - "))

### use not cln one, and add a row cbind rn and Rid ####
OGD_low_high_stat_gut <- OGD_low_high_stat_gut %>%
  mutate(rxnInfo = paste(rxn, rn, sep = " - "))

########------------------------------ Oral plot ----------------------------------#######
# Calculate the difference between Low and High mean values
OGD_diff <- OGD_low_high_stat_oral %>%
  mutate(Difference = low_mean - high_mean)

# Reorder rxnInfo based on the difference (from most enriched in Low to least)
OGD_diff <- OGD_diff %>%
  arrange(Difference) %>%
  mutate(rxnInfo = reorder(rxnInfo, Difference))

# Convert data to long format for easier plotting
OGD_long <- OGD_diff %>%
  select(rxnInfo, high_mean, low_mean) %>%
  pivot_longer(cols = c(high_mean, low_mean), names_to = "Condition", values_to = "Mean_Value")

# Ensure the x-axis is ordered according to OGD_diff$rxnInfo
OGD_long$rxnInfo <- factor(OGD_long$rxnInfo, levels = c(OGD_diff$rxnInfo))
levels(OGD_long$rxnInfo)

# Create a bar plot with log scale and ordered x-axis
# Create the bar plot
oral_plot <- ggplot(OGD_long, aes(x = rxnInfo, y = Mean_Value, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.85) + 
  scale_y_log10() +  # Use log scale due to large differences 
  labs(
    # title = "Comparison of High and Low Mean Values for Each rxn",
    x = "Reactobiome",
    y = "Mean Abundance Value (log scale)"
  ) + 
  scale_fill_manual(values = c("low_mean" = "#809FFF", "high_mean" = "#00008B")) +  # Custom colors 
  coord_flip() +  # Flip the coordinates to make the bars vertical 
  theme(
    legend.position = "top",           # Move the legend to the top 
    legend.title = element_blank(),    # Remove legend title 
    legend.text = element_text(size = 15), 
    axis.text.y = element_text(size = 15)  # Increase y-axis text size
  )

oral_plot # save oral_Rid_barplot_Low_vs_H 8*7 L


########------------------------------ Gut plot ----------------------------------#######
# Calculate the difference between Low and High mean values
OGD_diff_gut <- OGD_low_high_stat_gut %>%
  mutate(Difference = low_mean - high_mean)

# Reorder rxnInfo based on the difference (from most enriched in Low to least)
OGD_diff_gut <- OGD_diff_gut %>%
  arrange(Difference) %>%
  mutate(rxnInfo = reorder(rxnInfo, Difference))

# Convert data to long format for easier plotting
OGD_long_gut <- OGD_diff_gut %>%
  select(rxnInfo, high_mean, low_mean) %>%
  pivot_longer(cols = c(high_mean, low_mean), names_to = "Condition", values_to = "Mean_Value")

# Ensure the x-axis is ordered according to OGD_diff_gut$rxnInfo
OGD_long_gut$rxnInfo <- factor(OGD_long_gut$rxnInfo, levels = c(OGD_diff_gut$rxnInfo))
levels(OGD_long_gut$rxnInfo)

# Create a bar plot with log scale and ordered x-axis
# Create the bar plot
gut_plot <- ggplot(OGD_long_gut, aes(x = rxnInfo, y = Mean_Value, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.85) + 
  scale_y_log10() +  # Use log scale due to large differences 
  labs(
    # title = "Comparison of High and Low Mean Values for Each rxn",
    x = "Reactobiome",
    y = "Mean Abundance Value (log scale)"
  ) + 
  scale_fill_manual(values = c("low_mean" = "#809FFF", "high_mean" = "#00008B")) +  # Custom colors 
  coord_flip() +  # Flip the coordinates to make the bars vertical 
  theme(
    legend.position = "top",           # Move the legend to the top 
    legend.title = element_blank(),    # Remove legend title 
    legend.text = element_text(size = 15), 
    axis.text.y = element_text(size = 15)  # Increase y-axis text size
  )

gut_plot # save gut_Rid_barplot_Low_vs_H 9*11 P



