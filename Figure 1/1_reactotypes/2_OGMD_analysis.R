# load reative abundance reactobiome

reaction_relab_frommsp_gut_lc_and_h 
reaction_relab_frommsp_oral_lc_and_h

gut_rxty_pcoa_tab
oral_rxty_pcoa_tab

########################################################################################
############################# OGMD: Bray-Curtis Distance ###############################
########################################################################################
library(vegan)
library(ape)
library(phangorn)

rownames(rxbiome_rela_oral.t.m) = rxbiome_rela_oral.t.m[,1]

distance_matrix_gut <- vegdist(rxbiome_rela_gut.t.m[,-c(1:8)], method = "bray")
distance_matrix_gut.m = distance_matrix_gut %>% as.matrix() %>% as.data.frame()

distance_matrix_oral <- vegdist(rxbiome_rela_oral.t.m[,-c(1:8)], method = "bray")
distance_matrix_oral.m = distance_matrix_oral %>% as.matrix() %>% as.data.frame()

### merge metadata -- bray curtis distance ###
distance_matrix_oral.m$`SampleID` = rownames(distance_matrix_oral.m)
distance_matrix_gut.m$`SampleID` = rownames(distance_matrix_gut.m)

distance_matrix_oral.m.meta = merge(distance_matrix_oral.m,rxbiome_rela_oral.t.m[,1:8],by = 'SampleID', all.x=T)
distance_matrix_gut.m.meta = merge(distance_matrix_gut.m,rxbiome_rela_gut.t.m[,1:8],by = 'SampleID', all.x=T)

# PERMANOVA test
# Perform adonis test
adonis2(distance_matrix_oral ~ reactotype, data = distance_matrix_oral.m.meta, permutations = 999)
adonis2(distance_matrix_gut ~ reactotype, data = distance_matrix_gut.m.meta, permutations = 999)

adonis2(distance_matrix_oral ~ Severity_group, data = distance_matrix_oral.m.meta, permutations = 999)
adonis2(distance_matrix_gut ~ Severity_group, data = distance_matrix_gut.m.meta, permutations = 999)

adonis2(distance_matrix_oral ~ Group, data = distance_matrix_oral.m.meta, permutations = 999)
adonis2(distance_matrix_gut ~ Group, data = distance_matrix_gut.m.meta, permutations = 999)

adonis2(distance_matrix_oral ~ Gender, data = distance_matrix_oral.m.meta, permutations = 999)
adonis2(distance_matrix_gut ~ Gender, data = distance_matrix_gut.m.meta, permutations = 999)

adonis2(distance_matrix_oral ~ Age, data = distance_matrix_oral.m.meta, permutations = 999)
adonis2(distance_matrix_gut ~ Age, data = distance_matrix_gut.m.meta, permutations = 999)

####################
############################## bray distance oral gut #############################
# Step 1: Load necessary packages and data
library(vegan)

oral_matrix = rxbiome_rela_paired_oral_merged[,-c(1,111)]
gut_matrix = rxbiome_rela_paired_gut_merged[,-c(1,111)]
oral_matrix[is.na(oral_matrix)]=0
gut_matrix[is.na(gut_matrix)]=0

# Initialize an empty list to store the distances
paired_distances <- list()

# Loop through each sample and calculate the paired Bray distance
for (i in 1:109) {
  oral_sample <- oral_matrix[, i]
  gut_sample <- gut_matrix[, i]
  paired_distances[[i]] <- vegdist(rbind(oral_sample, gut_sample), method = "bray")[1]
}

# Create a dataframe with the paired distances and sample IDs
paired_distances_df <- data.frame(
  Paired_Distance = unlist(paired_distances),
  Oral_Sample_ID = sub("\\..*", "", colnames(oral_matrix)),
  Gut_Sample_ID = sub("\\..*", "", colnames(gut_matrix)))

paired_distances_df.m = merge(paired_distances_df, rxbiome_rela_gut.t.m[,1:8], by.x = "Gut_Sample_ID", by.y = "SampleID", all.x = TRUE)
paired_distances_df.m.oral = merge(paired_distances_df, rxbiome_rela_oral.t.m[,1:8], by.x = "Oral_Sample_ID", by.y = "SampleID", all.x = TRUE)

# ggplot
# severity group
ggplot(paired_distances_df.m, aes(x = Severity_group, y = Paired_Distance, fill = Severity_group)) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  labs(x = "Severity", y = "Oral-gut Rxbiome Bray distance", fill = "Severity group")+
  ggtitle("Boxplot of oral-gut distance by severity")+
  scale_fill_manual(values = severity_color, labels = c("Healthy","Mild Severity (3-8)", "Low Severity (9-14)", "Moderate Severity (15-24)", "High Severity (≥25)")) +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = Severity_group), label.y = 0.95)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(comparisons = list(
    # c("Healthy","Mild Severity"),
    # c("Healthy", "Low Severity"),
    # c("Healthy", "Moderate Severity"),
    c("Moderate Severity", "High Severity"),
    c("Healthy", "High Severity")
  ), label = "p.format") # p.format p.signif


### linear MELD ~ oral-gut_distance
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


