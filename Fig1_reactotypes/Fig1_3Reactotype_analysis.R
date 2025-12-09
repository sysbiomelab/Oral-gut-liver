# library(readr)

reactobiome_frommsp_gut_lc_and_h <- read.csv("~/Documents/.../reactobiome_frommsp_gut_lc_and_h_noFMT.csv",row.names = 1)
# 3608  169 # noFMT [1] 3608  146
reactobiome_frommsp_oral_lc_and_h <- read.csv("~/Documents/.../reactobiome_frommsp_oral_lc_and_h_noFMT.csv",row.names = 1)
# 2772  125 # noFMT [1] 2772  109

library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(parallel)

#### in create ####
# require(DirichletMultinomial)
# require(dplyr)
# require(parallel)
#### DMN clustering ###
# DMN for gut all
# set.seed(123)
# Fit_gutall <- mclapply(1:5, dmn, count=t(reactobiome_frommsp_gut_lc_and_h), verbose=TRUE)
# lplc_gutall <- sapply(Fit_gutall, laplace)
# best_gutall <- Fit_gutall[[which.min(lplc_gutall)]]


###### create finished dmn for gut based on t(reactobiome_frommsp_gut_lc_and_h)
load("~/Documents/.../dmm_gut.RData")

plot(lplc_gutall, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit")

p0 <- fitted(Fit_gutall[[1]], scale=TRUE) 
p3 <- fitted(best_gutall, scale=TRUE)
colnames(p3) <- paste("m", 1:3, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))
diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 40)
df_full <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff),3608)

heatmapdmn(t(reactobiome_frommsp_gut_lc_and_h), Fit_gutall[[1]], best_gutall, 40, lblwidth = 30)

heatmapdmn(t(reactobiome_frommsp_gut_lc_and_h), Fit_gutall[[1]], best_gutall, 20, lblwidth = 60)

# save table ### df, df_full
library(writexl)
df_save = df %>% as.data.frame()
df_save$Rxn = rownames(df_save)

df_full_save = df_full %>% as.data.frame()
df_full_save$Rxn = rownames(df_full_save)

################### contribution ################
diff_cluster = abs(p3 - as.vector(p0))
sum_diff_cluster = colSums(diff_cluster)
o_1 <- order(diff_cluster[,1], decreasing=TRUE)
o_2 <- order(diff_cluster[,2], decreasing=TRUE)
o_3 <- order(diff_cluster[,3], decreasing=TRUE)

cdiff_1 = cumsum(diff_cluster[o_1,1])
cdiff_2 = cumsum(diff_cluster[o_2,2])
cdiff_3 = cumsum(diff_cluster[o_3,3])

cdiff.Ratio_1 = cdiff_1 / sum_diff_cluster[1]
cdiff.Ratio_2 = cdiff_2 / sum_diff_cluster[2]
cdiff.Ratio_3 = cdiff_3 / sum_diff_cluster[3]

df_1 = cbind(Mean = p0[o_1], m1 = p3[o_1,1], diff_1 = diff_cluster[o_1,1],cdiff.Ratio_1)
df_1_head = head(cbind(Mean = p0[o_1], m1 = p3[o_1,1], diff_1 = diff_cluster[o_1,1],cdiff.Ratio_1),40)
rownames(df_1_head)[1:10]
# "R00764" "R02073" "R00206" "R10649" "R00019" "R10652" "R02017" "R02018" "R02019" "R02024"
# "R00764" "R02073" "R00206" "R10649" "R00019" "R10652" "R02017" "R02018" "R02019" "R02024"

df_2 = cbind(Mean = p0[o_2], m2 = p3[o_2,2], diff_2 = diff_cluster[o_2,2],cdiff.Ratio_2)
df_2_head = head(df_2,40)
rownames(df_2_head)[1:10]
# "R02017" "R02018" "R02019" "R02024" "R11170" "R11172" "R01528" "R10221" "R00835" "R02736"
# "R02017" "R02018" "R02019" "R02024" "R11170" "R11172" "R01528" "R10221" "R00835" "R02736"

df_3 = cbind(Mean = p0[o_3], m3 = p3[o_3,3], diff_3 = diff_cluster[o_3,3],cdiff.Ratio_3)
df_3_head = head(df_3,40)
rownames(df_3_head)[1:10]
# "R00764" "R02073" "R00206" "R00019" "R10649" "R02017" "R02018" "R02019" "R02024" "R00257"
# "R00764" "R02073" "R00206" "R00019" "R10649" "R02017" "R02018" "R02019" "R02024" "R00257"

rownames(df)[1:10]
# "R00764" "R02073" "R00206" "R02017" "R02018" "R02019" "R02024" "R10649" "R00019" "R10652"
# "R00764" "R02073" "R00206" "R02017" "R02018" "R02019" "R02024" "R10649" "R00019" "R10652"
#################################################################################


############ assign cluster #############
clusterAssigned = apply(Fit_gutall[[3]]@group, 1, function(x) which.max(x))

clusterAssigned.gut.rxtypes = clusterAssigned
clusterAssigned.gut.rxtypes = clusterAssigned.gut.rxtypes %>% as.data.frame()
colnames(clusterAssigned.gut.rxtypes)="reactotype"

##################################################################
################# load dmn_oral.RData ################
load("~/Documents/.../dmm_oral.RData")
# Fit_oralall <- mclapply(1:5, dmn, count=t(reactobiome_frommsp_oral_lc_and_h), verbose=TRUE)
# lplc_oralall <- sapply(Fit_oralall, laplace)
# best_oralall <- Fit_oralall[[which.min(lplc_oralall)]] # k=2
plot(lplc_oralall, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit") # saved 4*5
# save Model_fit_oral_reto 4*5
#############################

p0_o <- fitted(Fit_oralall[[1]], scale=TRUE) 
p3_o <- fitted(best_oralall, scale=TRUE)
colnames(p3_o) <- paste("m", 1:2, sep="")
(meandiff_o <- colSums(abs(p3_o - as.vector(p0_o))))
diff_o <- rowSums(abs(p3_o - as.vector(p0_o)))
o_o <- order(diff_o, decreasing=TRUE)
cdiff_o <- cumsum(diff_o[o_o]) / sum(diff_o)
df_o <- head(cbind(Mean=p0_o[o_o], p3_o[o_o,], diff=diff_o[o_o], cdiff_o), 40)
df_full_o <- cbind(Mean=p0_o[o_o], p3_o[o_o,], diff=diff_o[o_o], cdiff_o)

heatmapdmn(t(reactobiome_frommsp_oral_lc_and_h), Fit_oralall[[1]], best_oralall, 40, lblwidth = 60)
heatmapdmn(t(reactobiome_frommsp_oral_lc_and_h), Fit_oralall[[1]], best_oralall, 20, lblwidth = 60)
# save Heatmap_Rid_oral_reto 5*5

# save table ### df_o, df_full_o
library(writexl)
df_o_save = df_o %>% as.data.frame()
df_o_save$Rxn = rownames(df_o_save)
write_xlsx(df_o_save, "~/.../Overall_diff_oral_40.xlsx")

df_full_o_save = df_full_o %>% as.data.frame()
df_full_o_save$Rxn = rownames(df_full_o_save)
write_xlsx(df_full_o_save, "~/.../Overall_diff_oral_full.xlsx")

###########################################################################
diff_cluster_o = abs(p3_o - as.vector(p0_o))
sum_diff_cluster_o = colSums(diff_cluster_o)
o_1_o <- order(diff_cluster_o[,1], decreasing=TRUE)
o_2_o <- order(diff_cluster_o[,2], decreasing=TRUE)

cdiff_1_o = cumsum(diff_cluster_o[o_1_o,1])
cdiff_2_o = cumsum(diff_cluster_o[o_2_o,2])


cdiff.Ratio_1_o = cdiff_1_o / sum_diff_cluster_o[1]
cdiff.Ratio_2_o = cdiff_2_o / sum_diff_cluster_o[2]

df_1_o = cbind(Mean = p0_o[o_1_o],m1 = p3_o[o_1_o,1], diff_1 = diff_cluster_o[o_1_o,1],cdiff.Ratio_1_o)
df_1_head_o = head(df_1_o,40)
rownames(df_1_head_o)[1:10]
# "R00206" "R10649" "R03254" "R10652" "R09030" "R12500" "R01800" "R02302" "R11098" "R11099"
# noFMT: [1] "R00206" "R10649" "R03254" "R09030" "R10652" "R11098" "R11099" "R01175" "R04751" "R02661"

df_2_o = cbind(Mean = p0_o[o_2_o],m2 = p3_o[o_2_o,2], diff_2 = diff_cluster_o[o_2_o,2],cdiff.Ratio_2_o)
df_2_head_o = head(df_2_o,40)
rownames(df_2_head_o)[1:10]
# "R00206" "R10649" "R03254" "R01737" "R08878" "R12500" "R00257" "R02302" "R09030" "R10220"
# noFMT: [1] "R00206" "R10649" "R03254" "R01737" "R00257" "R08878" "R03470" "R03165" "R11098" "R11099"

rownames(df_o)[1:10]
# "R00206" "R10649" "R03254" "R12500" "R01737" "R08878" "R02302" "R09030" "R10652" "R10220"
# noFMT: [1] "R00206" "R10649" "R03254" "R08878" "R00257" "R03470" "R01737" "R11098" "R11099" "R09030"


########################################################################################
########### save files ###############
# save table ### df_1, df_2, df_3
library(writexl)

# gut
df_1_save = df_1 %>% as.data.frame()
df_1_save$Rxn = rownames(df_1_save)
write_xlsx(df_1_save, "~/Documents/.../diff_gut_R1.xlsx")

df_2_save = df_2 %>% as.data.frame()
df_2_save$Rxn = rownames(df_2_save)
write_xlsx(df_2_save, "~/Documents/.../diff_gut_R2.xlsx")

df_3_save = df_3 %>% as.data.frame()
df_3_save$Rxn = rownames(df_3_save)
write_xlsx(df_3_save, "~/Documents/.../diff_gut_R3.xlsx")

# oral
df_1_o_save = df_1_o %>% as.data.frame()
df_1_o_save$Rxn = rownames(df_1_o_save)
write_xlsx(df_1_o_save, "~/.../diff_oral_R1.xlsx")

df_2_o_save = df_2_o %>% as.data.frame()
df_2_o_save$Rxn = rownames(df_2_o_save)
write_xlsx(df_2_o_save, "~/.../diff_oral_R2.xlsx")


###############################
# load metadata for samples
library(readxl)
lcgutrxtotype_metadata <- read_excel("~/Documents/.../rxtotype_metadata_noFMT.xlsx", 
                                     sheet = "LC_gut")
loralcrxtotype_metadata <- read_excel("~/Documents/.../rxtotype_metadata_noFMT.xlsx", 
                                      sheet = "LC_oral")
allgutrxtotype_metadata <- read_excel("~/Documents/.../rxtotype_metadata_noFMT.xlsx", 
                                      sheet = "gut_all")
alloralrxtotype_metadata <- read_excel("~/Documents/.../rxtotype_metadata_noFMT.xlsx", 
                                       sheet = "oral_all")

library(dplyr)
## all
allgutrxtotype_metadata <- allgutrxtotype_metadata %>%
  mutate(Severity_group = case_when(
    MELD >= 3 & MELD <= 8 ~ "Mild Severity",
    MELD > 8 & MELD <= 14 ~ "Low Severity",
    MELD > 14 & MELD <= 24 ~ "Moderate Severity",
    MELD > 24 ~ "High Severity",
    TRUE ~ "Healthy"
  ))

## gut
lcgutrxtotype_metadata <- lcgutrxtotype_metadata %>%
  mutate(Severity_group = case_when(
    MELD >= 3 & MELD <= 8 ~ "Mild Severity",
    MELD > 8 & MELD <= 14 ~ "Low Severity",
    MELD > 14 & MELD <= 24 ~ "Moderate Severity",
    MELD > 24 ~ "High Severity",
    TRUE ~ NA_character_
  ))

# View the modified data frame
head(lcgutrxtotype_metadata)

## oral
loralcrxtotype_metadata <- loralcrxtotype_metadata %>%
  mutate(Severity_group = case_when(
    MELD >= 3 & MELD <= 8 ~ "Mild Severity",
    MELD > 8 & MELD <= 14 ~ "Low Severity",
    MELD > 14 & MELD <= 24 ~ "Moderate Severity",
    MELD > 4 ~ "High Severity",
    TRUE ~ NA_character_
  ))

head(loralcrxtotype_metadata)
##### PCoA #########
library(ape)
library(vegan)

gut.rxty.dist = vegdist(t(reactobiome_frommsp_gut_lc_and_h), method = "bray")
gut.rxty.pcoa = pcoa(gut.rxty.dist)
gut.rxty.pca = prcomp(t(reactobiome_frommsp_gut_lc_and_h))

order_gut.pcoa = rownames(gut.rxty.pcoa$vectors)
gut.cols.df = allgutrxtotype_metadata[match(order_gut.pcoa, allgutrxtotype_metadata$SampleID), ]

gut.cols.df$color <- ifelse(gut.cols.df$Severity_group == "Mild Severity", "skyblue",
                            ifelse(gut.cols.df$Severity_group == "Low Severity", "cyan",
                                   ifelse(gut.cols.df$Severity_group == "Moderate Severity", "blue",
                                          ifelse(gut.cols.df$Severity_group == "High Severity", "purple", "grey"))))


plot(gut.rxty.pcoa$vectors[,1:2], pch=16, col=gut.cols.df$color)
plot(gut.rxty.pca$x[,1:2], pch=16, col=gut.cols.df$color)


gut.rxty.pcoa.tab = data.frame(gut.rxty.pcoa$vectors[,1:2], 
                               class = rep("gray", dim(gut.rxty.pcoa$vectors)[1]))
gut.rxty.pcoa.tab = cbind(gut.rxty.pcoa.tab,gut.cols.df)
gut.rxty.pcoa.tab$Severity_group <- factor(gut.rxty.pcoa.tab$Severity_group, levels = c("Healthy","Mild Severity", "Low Severity", "Moderate Severity", "High Severity"))

clusterAssigned.gut.rxtypes$SampleID = rownames(clusterAssigned.gut.rxtypes)
gut.rxty.pcoa.tab = merge(gut.rxty.pcoa.tab,clusterAssigned.gut.rxtypes,by="SampleID",all.x=T)

gut.rxty.pcoa.tab$clusters = clusterAssigned.gut.rxtypes[match(gut.rxty.pcoa.tab$SampleID, rownames(clusterAssigned.gut.rxtypes)),1]

## class transfor
gut.rxty.pcoa.tab$Axis.1 <- as.numeric(gut.rxty.pcoa.tab$Axis.1)
gut.rxty.pcoa.tab$Axis.2 <- as.numeric(gut.rxty.pcoa.tab$Axis.2)

# Convert columns to factors
gut.rxty.pcoa.tab$clusters <- as.factor(gut.rxty.pcoa.tab$clusters)
gut.rxty.pcoa.tab$reactotype <- as.factor(gut.rxty.pcoa.tab$reactotype)
gut.rxty.pcoa.tab$Severity_group <- as.factor(gut.rxty.pcoa.tab$Severity_group)

# Check the structure of the data frame
str(gut.rxty.pcoa.tab)

library(ggplot2)
ggplot(gut.rxty.pcoa.tab, aes(x=Axis.1, y=Axis.2, color=Severity_group,shape=reactotype)) + 
  geom_point(size=1.5) + 
  #scale_color_manual(values = c("blue","cyan","gray","purple","skyblue")) + #stat_ellipse()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        #legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# gut reactitype PCoA plot
ggplot(gut.rxty.pcoa.tab, aes(x=Axis.1, y=Axis.2, color=reactotype,shape=Severity_group)) + 
  geom_point(size=1.5) + 
  scale_color_manual(values = reactotype_colors) + # stat_ellipse()+
  labs(title = "PCoA of Gut Reactobiome", x = "PCoA1", y = "PCoA2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        # legend.position = "top",
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))


split(gut.rxty.pcoa.tab$SampleID, gut.rxty.pcoa.tab$clusters)
lapply(split(gut.rxty.pcoa.tab$SampleID, gut.rxty.pcoa.tab$clusters), length)
# sample size of Severity_group
lapply(split(gut.rxty.pcoa.tab$SampleID, gut.rxty.pcoa.tab$Severity_group), length)
# sample size, gut samples
library(gplots)
par(mfrow = c(1, 1), mar = c(5, 6, 4, 2) + 0.1)
barplot2(table(gut.rxty.pcoa.tab$Severity_group), 
         col=c("gray","cyan","skyblue","blue","purple"),las=2,
         ylim=c(0,70),xlab="Severity Group", ylab="Count", border=NA,
         names.arg=c("Healthy","Mild Severity", "Low Severity", "Moderate Severity", "High Severity"))
axis(1, at=1:5, labels=c("Healthy","Mild Severity", "Low Severity", "Moderate Severity", "High Severity"), las=2, cex.axis=0.8)

ggplot(gut.rxty.pcoa.tab, aes(x = Severity_group)) +
  geom_bar(fill = c("gray", "cyan", "skyblue", "blue", "purple")) +
  geom_text(stat='count', aes(label=..count..), vjust=-0.3) +
  labs(x = "Severity Group", y = "Count") +
  ylim(0, 70)+
  ggtitle("Gut") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
# save Group_count_gut.pdf 3*4 P

# Create the grouped bar plot to visualize the distribution of severity groups across clusters
ggplot(gut.rxty.pcoa.tab, aes(x = as.factor(reactotype), fill = Severity_group)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = ..count.., group = Severity_group), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(x = "Reactotype", y = "Count", fill = "Severity Group") +
  ylim(0, 72)+
  ggtitle("Distribution of Severity Groups Across Gut Reactotypes") +
  scale_fill_manual(values = c("gray", "cyan", "skyblue", "blue", "purple"),
                    labels = c("Healthy", "Mild Severity", "Low Severity", "Moderate Severity", "High Severity"))+
  theme_minimal() # reactotype_group_count.pdf 5.5*4 L

## sig test
# Calculate the mean and standard deviation for age, MELD, and BMI for each reactotype
gut.rxty.pcoa.tab$Age = gut.rxty.pcoa.tab$Age %>% as.numeric()
gut.rxty.pcoa.tab$BMI = gut.rxty.pcoa.tab$BMI %>% as.numeric()
gut.rxty.pcoa.tab$MELD = gut.rxty.pcoa.tab$MELD %>% as.numeric()

# table for summerise the Mean and SD for Age BMI MELD in each reactotype in GUT
summary_df_gut_lc <- aggregate( . ~ reactotype, 
                                data = gut.rxty.pcoa.tab[-which(gut.rxty.pcoa.tab$Severity_group=="Healthy"), c("Age", "MELD", "BMI","reactotype")], 
                                function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE)))
summary_df_gut_h <- aggregate( . ~ reactotype, 
                               data = gut.rxty.pcoa.tab[which(gut.rxty.pcoa.tab$Severity_group=="Healthy"), c("Age", "BMI","reactotype")], 
                               function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE)))

# Create box plots for age, MELD, and BMI for each reactotype
library(ggpubr)
plot_age <- ggplot(gut.rxty.pcoa.tab, aes(x = reactotype, y = Age, fill = reactotype)) +
  geom_boxplot() +
  labs(x = "Reactotype", y = "Age") +
  ggtitle("Box Plot for Age by Reactotype")

# Manually set colors for reactotype levels
reactotype_colors <- c("#B0E0E6", "#5F9EA0", "#1F78B4") # 
# c("#ADD8E6", "#4169E1", "#000080") # c("#66C2A5", "#FC8D62", "#8DA0CB")  # Example of colors with a gradient from low to high

# Age
# Create box plots for Age by Reactotype
ggplot(gut.rxty.pcoa.tab[-which(gut.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = reactotype, y = Age, fill = reactotype)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list( c("1", "3")), label = "p.signif", hide.ns = TRUE) +
  labs(x = "Gut Reactotype", y = "Age") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Box Plot for Age by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)

# Create box plots for Age by Reactotype
ggplot(gut.rxty.pcoa.tab[which(gut.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = reactotype, y = Age, fill = reactotype)) +
  geom_boxplot() +
  #stat_compare_means(comparisons = list(), label = "p.signif", hide.ns = TRUE) +
  labs(x = "Gut Reactotype", y = "Age") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Box Plot for Age by Reactotype of Healthy") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)

# Create box plots for Age by Reactotype
ggplot(gut.rxty.pcoa.tab[,], aes(x = reactotype, y = Age, fill = Group)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("H","LC")), label = "p.signif", hide.ns = TRUE) +
  labs(x = "Gut Reactotype", y = "Age") +
  scale_fill_manual(values = c("grey","#8DA0CB")) +
  #facet_wrap(~ reactotype, scales = "free", ncol = 1) +
  ggtitle("Box Plot for Age by Reactotype HC, LC") +
  theme_minimal()

# Convert the 'reactotype' column to a factor
gut.rxty.pcoa.tab$reactotype <- as.factor(gut.rxty.pcoa.tab$reactotype)

# Verify the data type
class(gut.rxty.pcoa.tab$reactotype)

ggplot(gut.rxty.pcoa.tab, aes(x=reactotype, y=Age)) + 
  geom_boxplot()

## keep digits for BMI
gut.rxty.pcoa.tab$BMI <- round(gut.rxty.pcoa.tab$BMI, 2)

## keep digits for BMI
gut.rxty.pcoa.tab$MELD <- round(gut.rxty.pcoa.tab$MELD)

# MELD
# Create box plots for Age by Reactotype
ggplot(gut.rxty.pcoa.tab[-which(gut.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = reactotype, y = MELD, fill = reactotype)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")), label = "p.signif") +
  labs(x = "Gut Reactotype", y = "MELD") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Box Plot for MELD by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)

# BMI
# Create box plots for Age by Reactotype
ggplot(gut.rxty.pcoa.tab[-which(gut.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = reactotype, y = BMI, fill = reactotype)) +
  geom_boxplot() +
  # stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")), label = "p.signif") +
  labs(x = "Gut Reactotype", y = "BMI") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Box Plot for BMI by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)
# BMI HC
ggplot(gut.rxty.pcoa.tab[which(gut.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = reactotype, y = BMI, fill = reactotype)) +
  geom_boxplot() +
  #stat_compare_means(comparisons = list(c("1", "2")), label = "p.signif") +
  labs(x = "Gut Reactotype", y = "BMI") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Box Plot for BMI by Reactotype of HC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)

# Plot the stacked bar plot
ggplot(gut.rxty.pcoa.tab, aes(x = Severity_group, fill = as.factor(reactotype))) +
  geom_bar(position = "fill") +
  labs(x = "Severity Group", y = "Count", fill = "Reactotype") +
  ggtitle("Reactotypes in Each Severity Group") +
  scale_fill_manual(values = reactotype_colors) +  # Customize colors as needed
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top") # Reactotypes in Each Severity Group gut proportion 3*4

ggplot(gut.rxty.pcoa.tab, aes(x = Severity_group, fill = as.factor(reactotype))) +
  geom_bar(position = "stack") +
  labs(x = "Severity Group", y = "Count", fill = "Reactotype") +
  ggtitle("Reactotypes in Each Severity Group") +
  scale_fill_manual(values = reactotype_colors) +  # Customize colors as needed
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top") # Reactotypes in Each Severity Group gut count 3*4

#
# Plot the stacked bar plot
ggplot(gut.rxty.pcoa.tab, aes(x = reactotype, fill = as.factor(Severity_group))) +
  geom_bar(position = "fill") +
  labs(x = "Severity Group", y = "Count", fill = "Reactotype") +
  ggtitle("Relative Proportions Gut Reactotype") +
  scale_fill_manual(values = c("gray", "cyan", "skyblue", "blue", "purple")) +  # Customize colors as needed
  theme_minimal()
# +theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top") # Reactotype_gut_composition_relative.pdf 4*4

ggplot(gut.rxty.pcoa.tab, aes(x = reactotype , fill = as.factor(Severity_group))) +
  geom_bar(position = "stack") +
  labs(x = "Severity Group", y = "Count", fill = "Reactotype") +
  ggtitle("Stack Barplot Gut reactotype") +
  scale_fill_manual(values = c("gray", "cyan", "skyblue", "blue", "purple")) +  # Customize colors as needed
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()
+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top") 



# raincould plots

# MELD
# Create box plots for Age by Reactotype
library(tidyverse)
library(palmerpenguins)
library(ggdist)
packageVersion("ggdist")

ggplot(gut.rxty.pcoa.tab[-which(gut.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = reactotype, y = MELD, fill = reactotype)) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")), label = "p.format") + # p.format p.signif
  labs(x = "Gut Reactotype", y = "MELD") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Raincloud for MELD by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = reactotype), label.y = 53)
+geom_jitter(shape=1, position = position_jitter(0.2))
# Raincloud_for_MELD_by_Reactotype_of_LC_gut 4*4

# Age 
ggplot(gut.rxty.pcoa.tab[,], aes(x = reactotype, y = Age, fill = reactotype)) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_compare_means(comparisons = list( c("1", "3")), label = "p.format", hide.ns = TRUE) + #p.signif
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  labs(x = "Gut Reactotype", y = "Age") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Raincloud for Age by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = reactotype), label.y = 90)
# +geom_jitter(shape=1, position = position_jitter(0.2))
# Raincloud_for_Age_by_Reactotype_of_LC_gut 4*4

# BMI 
ggplot(gut.rxty.pcoa.tab[,], aes(x = reactotype, y = BMI, fill = reactotype)) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  labs(x = "Gut Reactotype", y = "BMI") +
  scale_fill_manual(values = reactotype_colors) +
  ggtitle("Raincloud for BMI by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = reactotype))
# +geom_jitter(shape=1, position = position_jitter(0.2))
# Raincloud_for_BMI_by_Reactotype_of_LC_gut 4*4






##### oral
##### PCoA #########
library(ape)
oral.rxty.dist = vegdist(t(reactobiome_frommsp_oral_lc_and_h), method = "bray")
oral.rxty.pcoa = pcoa(oral.rxty.dist)
oral.rxty.pca = prcomp(t(reactobiome_frommsp_oral_lc_and_h))

order_oral.pcoa = rownames(oral.rxty.pcoa$vectors)

## oral
alloralrxtotype_metadata <- alloralrxtotype_metadata %>%
  mutate(Severity_group = case_when(
    MELD >= 3 & MELD <= 10 ~ "Mild Severity",
    MELD > 10 & MELD <= 18 ~ "Low Severity",
    MELD > 18 & MELD <= 28 ~ "Moderate Severity",
    MELD > 28 ~ "High Severity",
    TRUE ~ "Healthy"
  ))

head(alloralrxtotype_metadata)

oral.cols.df = alloralrxtotype_metadata[match(order_oral.pcoa, alloralrxtotype_metadata$SampleID), ]



oral.cols.df$color <- ifelse(oral.cols.df$Severity_group == "Mild Severity", "skyblue",
                             ifelse(oral.cols.df$Severity_group == "Low Severity", "cyan",
                                    ifelse(oral.cols.df$Severity_group == "Moderate Severity", "blue",
                                           ifelse(oral.cols.df$Severity_group == "High Severity", "purple", "grey"))))


plot(oral.rxty.pcoa$vectors[,1:2], pch=16, col=oral.cols.df$color)
plot(oral.rxty.pca$x[,1:2], pch=16, col=oral.cols.df$color)


oral.rxty.pcoa.tab = data.frame(oral.rxty.pcoa$vectors[,1:2])
oral.rxty.pcoa.tab = cbind(oral.rxty.pcoa.tab,oral.cols.df)
oral.rxty.pcoa.tab$Severity_group <- factor(oral.rxty.pcoa.tab$Severity_group, levels = c("Healthy","Mild Severity", "Low Severity", "Moderate Severity", "High Severity"))

clusterAssigned.oral.rxtypes$SampleID = rownames(clusterAssigned.oral.rxtypes)
oral.rxty.pcoa.tab = merge(oral.rxty.pcoa.tab,clusterAssigned.oral.rxtypes,by="SampleID",all.x=T)

oral.rxty.pcoa.tab$clusters = clusterAssigned.oral.rxtypes[match(oral.rxty.pcoa.tab$SampleID, rownames(clusterAssigned.oral.rxtypes)),1]

oral.rxty.pcoa.tab$reactotype = as.factor(oral.rxty.pcoa.tab$reactotype)
oral.rxty.pcoa.tab$clusters = as.factor(oral.rxty.pcoa.tab$clusters)

##
##  BMI
oral.rxty.pcoa.tab$BMI = oral.rxty.pcoa.tab$BMI %>% as.numeric()
oral.rxty.pcoa.tab$BMI <- round(oral.rxty.pcoa.tab$BMI, 2)

##  MELD
oral.rxty.pcoa.tab$MELD = oral.rxty.pcoa.tab$MELD %>% as.numeric()
oral.rxty.pcoa.tab$MELD <- round(oral.rxty.pcoa.tab$MELD)

# numeric Age
oral.rxty.pcoa.tab$Age = oral.rxty.pcoa.tab$Age %>% as.numeric()

#

head(oral.rxty.pcoa.tab)

##### oral plotting ######
# PCoA  --- Oral reactotype PCoA plot
ggplot(oral.rxty.pcoa.tab, aes(x=Axis.1, y=Axis.2, color=factor(reactotype),shape=Severity_group)) + 
  geom_point(size=1.5) + 
  scale_color_manual(values = c("blue","red","gold")) + # stat_ellipse()+
  labs(title = "PCoA of Gut Reactobiome", x = "PCoA1", y = "PCoA2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        # legend.position = "top",
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))


# Bar plot ---- type - severity - Bar plot percentage
split(oral.rxty.pcoa.tab$SampleID, oral.rxty.pcoa.tab$reactotype)
lapply(split(oral.rxty.pcoa.tab$SampleID, oral.rxty.pcoa.tab$reactotype), length) # 1--108; 2--16 
# sample size of Severity_group
lapply(split(oral.rxty.pcoa.tab$SampleID, oral.rxty.pcoa.tab$Severity_group), length)
# $Healthy[1] 56 $`Mild Severity`[1] 16 $`Low Severity`[1] 36  $`Moderate Severity`[1] 13  $`High Severity`[1] 3
# sample size of oral samples
library(gplots)
par(mfrow = c(1, 1), mar = c(5, 6, 4, 2) + 0.1)

ggplot(oral.rxty.pcoa.tab, aes(x = Severity_group)) +
  geom_bar(fill = c("gray", "cyan", "skyblue", "blue", "purple")) +
  geom_text(stat='count', aes(label=..count..), vjust=-0.3) +
  labs(x = "Severity Group", y = "Count") +
  ylim(0, 70)+
  ggtitle("Oral histogram") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Group_count_oral.pdf 3*4 P

# Create the grouped bar plot to visualize the distribution of severity groups across clusters
ggplot(oral.rxty.pcoa.tab, aes(x = as.factor(reactotype), fill = Severity_group)) +
  geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = ..count.., group = Severity_group), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(x = "Reactotype", y = "Count", fill = "Severity Group") +
  ylim(0, 72)+
  ggtitle("Distribution of Severity Groups Across Reactotypes Oral") +
  scale_fill_manual(values = c("gray", "cyan", "skyblue", "blue", "purple"),
                    labels = c("Healthy", "Mild Severity", "Low Severity", "Moderate Severity", "High Severity"))+
  theme_minimal() # Distribution of Severity Groups Across Reactotypes Oral.pdf 5*4 L

# Plot the stacked bar plot
ggplot(oral.rxty.pcoa.tab, aes(x = reactotype, fill = as.factor(Severity_group))) +
  geom_bar(position = "fill") +
  labs(x = "Severity Group", y = "Count", fill = "Reactotype") +
  ggtitle("Relative Proportions Oral Reactotype") +
  scale_fill_manual(values = c("gray", "cyan", "skyblue", "blue", "purple")) +  # Customize colors as needed
  theme_minimal()
# +theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top") # Reactotype_oral_composition_relative.pdf 3.5*4 L

ggplot(oral.rxty.pcoa.tab, aes(x = reactotype , fill = as.factor(Severity_group))) +
  geom_bar(position = "stack") +
  labs(x = "Severity Group", y = "Count", fill = "Reactotype") +
  ggtitle("Stack Barplot Oral reactotype") +
  scale_fill_manual(values = c("gray", "cyan", "skyblue", "blue", "purple")) +  # Customize colors as needed
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()
+theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top") # Stack Barplot Oral reactotype 3.5*4

#  MELD oral
ggplot(oral.rxty.pcoa.tab[-which(oral.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = factor(reactotype), y = MELD, fill = factor(reactotype))) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  stat_compare_means(comparisons = list(c("1", "2")), label = "p.format") + # p.format p.signif
  labs(x = "Oral Reactotype", y = "MELD") +
  scale_fill_manual(values = reactotype_colors[1:2]) +
  ggtitle("Raincloud for MELD by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = reactotype), label.y = 34)
+geom_jitter(shape=1, position = position_jitter(0.2))
#save Raincloud_for_MELD_by_Reactotype_of_LC_oral 4*4

#  BMI oral
ggplot(oral.rxty.pcoa.tab[-which(oral.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = factor(reactotype), y = BMI, fill = factor(reactotype))) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  # stat_compare_means(comparisons = list(c("1", "2")), label = "p.signif") +
  labs(x = "Oral Reactotype", y = "BMI") +
  scale_fill_manual(values = reactotype_colors[1:2]) +
  ggtitle("Raincloud for BMI by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = reactotype), label.y = 55)
+geom_jitter(shape=1, position = position_jitter(0.2))
#save Raincloud_for_BMI_by_Reactotype_of_LC_oral 4*4

#  Age oral
ggplot(oral.rxty.pcoa.tab[-which(oral.rxty.pcoa.tab$Severity_group=="Healthy"),], aes(x = factor(reactotype), y = Age, fill = factor(reactotype))) +
  stat_slab(aes(thickness = stat(pdf*n)),
            scale = 0.7) +
  stat_dotsinterval(side = "bottom", 
                    scale = 0.7, 
                    slab_size = NA)+
  # stat_compare_means(comparisons = list(c("1", "2")), label = "p.signif") +
  labs(x = "Oral Reactotype", y = "Age") +
  scale_fill_manual(values = reactotype_colors[1:2]) +
  ggtitle("Raincloud for Age by Reactotype of LC") +
  theme_minimal()+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size=4)+
  stat_compare_means(aes(group = reactotype), label.y = 76)
+geom_jitter(shape=1, position = position_jitter(0.2))
#save Raincloud_for_Age_by_Reactotype_of_LC_oral 4*4

library(writexl)
write_xlsx(oral.rxty.pcoa.tab, "~/.../oral.rxty.pcoa.tab.xlsx")
write_xlsx(gut.rxty.pcoa.tab, "~/.../gut.rxty.pcoa.tab.xlsx")

