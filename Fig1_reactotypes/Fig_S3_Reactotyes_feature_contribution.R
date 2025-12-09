################################################
######### reactotype top contribution #####
################################################

###### when got finished dmn for gut based on t(reactobiome_frommsp_gut_lc_and_h)
# load("/Users/jinyi/Documents/.../dmn.RData")

reactobiome_frommsp_gut_lc_and_h <- read.csv("~/Documents/.../reactobiome_frommsp_gut_lc_and_h_noFMT.csv",row.names = 1)
# 3608  169 # noFMT [1] 3608  146
# reactobiome_frommsp_oral_lc_and_h <- read.csv("~/Documents/.../reactobiome_frommsp_oral_lc_and_h_noFMT.csv",row.names = 1)
# # 2772  125 # noFMT [1] 2772  109


plot(lplc_gutall, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit") # save 3.5*4

p0 <- fitted(Fit_gutall[[1]], scale=TRUE) 
p3 <- fitted(best_gutall, scale=TRUE)
colnames(p3) <- paste("m", 1:3, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))
diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 40)
df_full <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff),3608)

heatmapdmn(t(reactobiome_frommsp_gut_lc_and_h), Fit_gutall[[1]], best_gutall, 40, lblwidth = 60)

###################-----------------------------------------------------------###################

library(reshape2)
#Contribution of each taxon to each component
colr <- c("#B0E0E6", "#5F9EA0", "#1F78B4") # 

# for (k in seq(ncol(fitted(best_gutall)))) {
#   d <- melt(fitted(best_gutall))
#   colnames(d) <- c("OTU", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#     # Arrange OTUs by assignment strength
#     arrange(value) %>%
#     mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
#     # Only show the most important drivers
#     filter(abs(value) > quantile(abs(value), 0.95))     
#   #plot to subfigures in pdf file
#   
#   p <- ggplot(d, aes(x = OTU, y = value, fill="a")) +
#     geom_bar( stat = "identity") +
#     coord_flip() +
#     theme_minimal()+
#     scale_fill_manual(values= colr[k])+
#     theme(legend.position="none")

#-----------------------------------------------------------#

# k = 1
d <- melt(fitted(best_gutall))
colnames(d) <- c("Rxn", "cluster", "value")
d <- subset(d, cluster == 1) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(Rxn = factor(Rxn, levels = unique(Rxn))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.99))

#plot to subfigures in pdf file

ggplot(d, aes(x = Rxn, y = value, fill="a")) +
  geom_bar( stat = "identity") +
  coord_flip() +
  theme_minimal()+
  scale_fill_manual(values= colr[1])+
  theme(legend.position="none")

# add Rxn annotation #
library(readr)
pathway_reaction_list <- read_csv("~/Documents/github_reactobiome/reactobiome/rytoPJ_231013/RF_231128/pathway_reaction_list.csv")
Rxn_anno = as.data.frame(pathway_reaction_list)
colnames(Rxn_anno) = c("rn","Rxn","name")

d_anno = merge(d,Rxn_anno,by="Rxn",all.x=T)

library(dplyr)
d_anno <- d_anno %>%
  mutate(Anno = paste(Rxn, rn, sep = " - "))
d_anno_1 = d_anno
#plot to subfigures in pdf file
ggplot(d_anno_1, aes(x = Anno, y = value, fill="a")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=colr[1]) +
  # theme(legend.position="none") +
  theme(
    legend.position = "none",  # 移除图例
    axis.text.y = element_text(hjust = 0)  # 将 y 轴的文本标签左对齐
  ) +
  ggtitle("Gut Reactotype 1")  # Add the title "xxx" to the plot

# save Gut_reto_1_Barplot 6*8 P
#-----------------------------------------------------------#

# k = 2
d <- melt(fitted(best_gutall))
colnames(d) <- c("Rxn", "cluster", "value")
d <- subset(d, cluster == 2) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(Rxn = factor(Rxn, levels = unique(Rxn))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.99))

d_anno = merge(d,Rxn_anno,by="Rxn",all.x=T)

library(dplyr)
d_anno <- d_anno %>%
  mutate(Anno = paste(Rxn, rn, sep = " - "))

d_anno_2 = d_anno
#plot to subfigures in pdf file
ggplot(d_anno_2, aes(x = Anno, y = value, fill="a")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=colr[2]) +
  # theme(legend.position="none") +
  theme(
    legend.position = "none",  # 移除图例
    axis.text.y = element_text(hjust = 0)  # 将 y 轴的文本标签左对齐
  ) +
  ggtitle("Gut Reactotype 2")  # Add the title "xxx" to the plot

# save Gut_reto_2_Barplot 6*8 P
#-----------------------------------------------------------#

# k = 3
d <- melt(fitted(best_gutall))
colnames(d) <- c("Rxn", "cluster", "value")
d <- subset(d, cluster == 3) %>%
  # Arrange OTUs by assignment strength
  arrange(value) %>%
  mutate(Rxn = factor(Rxn, levels = unique(Rxn))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.99))

d_anno = merge(d,Rxn_anno,by="Rxn",all.x=T)

library(dplyr)
d_anno <- d_anno %>%
  mutate(Anno = paste(Rxn, rn, sep = " - "))

d_anno_3 = d_anno
#plot to subfigures in pdf file
ggplot(d_anno_3, aes(x = Anno, y = value, fill="a")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=colr[3]) +
  # theme(legend.position="none") +
  theme(
    legend.position = "none",  # 移除图例
    axis.text.y = element_text(hjust = 0)  # 将 y 轴的文本标签左对齐
  ) +
  ggtitle("Gut Reactotype 3")  

# save Gut_reto_3_Barplot 6*8 P

##################################################################
# Venn
# Assuming d_anno_1, d_anno_2, and d_anno_3 are your data frames and Rxn is the column of interest

setwd("~/Documents/Reactotype_2024_new/Reacto_contributor_240902")
# Load the necessary library
library(VennDiagram)

# Define the three lists
list1 <- d_anno_1$Rxn
list2 <- d_anno_2$Rxn
list3 <- d_anno_3$Rxn

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(
    "List 1" = list1,
    "List 2" = list2,
    "List 3" = list3
  ),
  category.names = c("G1", "G2", "G3"),
  filename = NULL,  # Use NULL to avoid saving to file automatically
  output = TRUE
)

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)

#
# Create the Venn diagram
venn.plot2 <- venn.diagram(
  x = list(
    "List 1" = list1,
    "List 2" = list2,
    "List 3" = list3
  ),
  category.names = c("G1", "G2", "G3"),
  filename = NULL,  # Use NULL to avoid saving to file automatically
  output = TRUE,
  fill = colr,  # Specify the colors for the Venn diagram
  alpha = 0.5,  # Set transparency of the circles
  cat.col = colr,  # Color of category names
  cat.cex = 1.5,  # Size of category names
  # cat.pos = c(-20, 20, 0),  # Position of category names
  cat.dist = c(0.05, 0.05, 0.05),  # Distance of category names from circles
  cex = 1.5,  # Size of numbers inside the Venn diagram
  lwd = 2  # Line width of the circles
) 

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot2) # save Venn_Reacotype_gut 3*3

#########
# Select the elements in G1 (list1) that are not in G2 (list2) or G3 (list3)
G1_unique1 <- setdiff(list1, union(list2, list3))
G1_unique1
# [1] "R05578" "R03665" "R01067" "R01641" "R01830" "R06590" "R03659" "R04773" "R00248" "R04112"

G1_unique2 <- setdiff(list2, union(list1, list3))
G1_unique2
# [1] "R00899" "R04951"

G1_unique3 <- setdiff(list3, union(list2, list1))
G1_unique3
# [1] "R04710" "R01150" "R00382" "R00480" "R01230" "R01231" "R08244" "R03664" "R00660"

#####------------------------------------------------###
ggplot(d_anno_1[which(d_anno_1$Rxn %in% G1_unique1),], aes(x = Anno, y = value, fill="a")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=colr[1]) +
  # theme(legend.position="none") +
  theme(
    legend.position = "none",  # 移除图例
    axis.text.y = element_text(hjust = 0)  # 将 y 轴的文本标签左对齐
  ) +
  ggtitle("Gut Reactotype 1")  # Add the title "xxx" to the plot


ggplot(d_anno_2[which(d_anno_2$Rxn %in% G1_unique2),], aes(x = Anno, y = value, fill="a")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=colr[2]) +
  # theme(legend.position="none") +
  theme(
    legend.position = "none",  # 移除图例
    axis.text.y = element_text(hjust = 0)  # 将 y 轴的文本标签左对齐
  ) +
  ggtitle("Gut Reactotype 2")  # Add the title "xxx" to the plot


ggplot(d_anno_3[which(d_anno_3$Rxn %in% G1_unique3),], aes(x = Anno, y = value, fill="a")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=colr[3]) +
  # theme(legend.position="none") +
  theme(
    legend.position = "none",  # 移除图例
    axis.text.y = element_text(hjust = 0)  # 将 y 轴的文本标签左对齐
  ) +
  ggtitle("Gut Reactotype 3")  # Add the title "xxx" to the plot




