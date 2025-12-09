library(readxl)
gut_all_samples_16MSP <- read_excel("~/.../gut_all_samples_16MSP.xlsx", 
                                    sheet = "MSP_16")


# Load necessary libraries
library(ggplot2)
library(igraph)
library(dplyr)

# 1. Load your data (example CSV file where rows are MSPs and columns are samples)
# Assuming your data is in a CSV file where first row and column are labels
abundance_data <- gut_all_samples_16MSP[,-1] %>% as.data.frame()
rownames(abundance_data) = gut_all_samples_16MSP$Species
abundance_data_t = t(abundance_data)
# 2. Calculate the correlation matrix between MSPs (species) across samples
# You can use Pearson or Spearman correlation depending on your data
cor_matrix_spearman <- cor(abundance_data.t, method = "spearman")  # or method = "spearman","pearson"
cor_matrix_pearson <- cor(abundance_data.t, method = "pearson")  # or method = "spearman","pearson"

# Initialize empty matrices to store correlation coefficients and p-values
cor_matrix_spearman <- matrix(NA, nrow = ncol(abundance_data_t), ncol = ncol(abundance_data_t))
cor_matrix_pearson <- matrix(NA, nrow = ncol(abundance_data_t), ncol = ncol(abundance_data_t))
p_values_spearman <- matrix(NA, nrow = ncol(abundance_data_t), ncol = ncol(abundance_data_t))
p_values_pearson <- matrix(NA, nrow = ncol(abundance_data_t), ncol = ncol(abundance_data_t))

# Compute the correlation and p-values for each pair of MSPs
for (i in 1:(ncol(abundance_data_t) - 1)) {
  for (j in (i + 1):ncol(abundance_data_t)) {
    # Spearman correlation and p-value
    cor_test_spearman <- cor.test(abundance_data_t[, i], abundance_data_t[, j], method = "spearman")
    cor_matrix_spearman[i, j] <- cor_test_spearman$estimate
    p_values_spearman[i, j] <- cor_test_spearman$p.value
    
    # Pearson correlation and p-value
    cor_test_pearson <- cor.test(abundance_data_t[, i], abundance_data_t[, j], method = "pearson")
    cor_matrix_pearson[i, j] <- cor_test_pearson$estimate
    p_values_pearson[i, j] <- cor_test_pearson$p.value
  }
}

# Convert the matrices to data frames for easier viewing
cor_matrix_spearman_df <- as.data.frame(cor_matrix_spearman)
cor_matrix_pearson_df <- as.data.frame(cor_matrix_pearson)
p_values_spearman_df <- as.data.frame(p_values_spearman)
p_values_pearson_df <- as.data.frame(p_values_pearson)


##### make table for netwrok ####
# Initialize an empty data frame for the edge list
edge_list <- data.frame(Source = character(),
                        Target = character(),
                        Correlation = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)

# Loop through each pair of MSPs to generate the edge list with correlation and p-value
for (i in 1:(ncol(abundance_data_t) - 1)) {
  for (j in (i + 1):ncol(abundance_data_t)) {
    
    # Get the correlation and p-value for Spearman
    cor_value_spearman <- cor_matrix_spearman[i, j]
    p_value_spearman <- p_values_spearman[i, j]
    
    # Add to the edge list
    edge_list <- rbind(edge_list, data.frame(Source = colnames(abundance_data_t)[i],
                                             Target = colnames(abundance_data_t)[j],
                                             Correlation = cor_value_spearman,
                                             p_value = p_value_spearman))
    
    # # Optionally, do the same for Pearson correlation if desired
    # cor_value_pearson <- cor_matrix_pearson[i, j]
    # p_value_pearson <- p_values_pearson[i, j]
    # 
    # edge_list <- rbind(edge_list, data.frame(Source = colnames(abundance_data_t)[i],
    #                                          Target = colnames(abundance_data_t)[j],
    #                                          Correlation = cor_value_pearson,
    #                                          p_value = p_value_pearson))
  }
}

# select the p value <= 0.05, filter the correlation
edge_list_sig = edge_list[which(edge_list$p_value <= 0.05),]

# add 16 Bugs name 
library(readxl)
MSP_16_info <- read_excel("~/.../gut_in_oral_overlap_no_H.xlsx", 
                          sheet = "co_enriched")

# # Merge the two dataframes based on the matching columns
edge_list_sig_MSPname <- merge(edge_list_sig, MSP_16_info[, c("Species", "species.x")],
                               by.x = "Source", by.y = "Species", all.x = TRUE)

colnames(edge_list_sig_MSPname)[5] = "Source_name"

edge_list_sig_MSPname <- merge(edge_list_sig_MSPname, MSP_16_info[, c("Species", "species.x")],
                               by.x = "Target", by.y = "Species", all.x = TRUE)

colnames(edge_list_sig_MSPname)[6] = "Target_name"

setwd("~/Documents/Reactotype_2024_new/ReDo_16MSP_network")
# Save the edge list to a CSV file
library(writexl)
write_xlsx(edge_list_sig_MSPname, "msp_coabundance_network.xlsx")


library(readxl)
Flux <- read_excel("~/.../gut_in_oral_overlap_no_H.xlsx", 
                   sheet = "FBA_FVA")
colnames(edge_list_sig_MSPname)
colnames(Flux)

edge_list_sig_MSPname_Flux <- merge(edge_list_sig_MSPname, Flux[, c("FVA_maxFlux: Ex_NH3", "FBA: Ex_NH3","species")],
                                    by.x = "Source_name", by.y = "species", all.x = TRUE)

colnames(edge_list_sig_MSPname_Flux)[7:8] = c("Source_FVA_NH3","Source_FBA_NH3")

# Save the edge list to a CSV file
library(writexl)
write_xlsx(edge_list_sig_MSPname_Flux, "msp_coabundance_network_flux.xlsx")





# 3. Optionally, filter out correlations that are not significant
# Example: keep correlations with an absolute value greater than 0.7
threshold <- 0.7
cor_matrix_filtered <- cor_matrix[abs(cor_matrix) > threshold]

# 4. Convert the correlation matrix into an edge list (for Cytoscape)
# Create a data frame with edges
edges <- data.frame()
for (i in 1:(nrow(cor_matrix_filtered)-1)) {
  for (j in (i+1):ncol(cor_matrix_filtered)) {
    if (!is.na(cor_matrix_filtered[i,j]) && abs(cor_matrix_filtered[i,j]) > threshold) {
      edges <- rbind(edges, data.frame(Source = rownames(cor_matrix_filtered)[i],
                                       Target = colnames(cor_matrix_filtered)[j],
                                       Weight = cor_matrix_filtered[i,j]))
    }
  }
}

# 5. Save the edge list as a CSV file for Cytoscape
write.csv(edges, "msp_coabundance_edges.csv", row.names = FALSE)

# Optional: visualize the network using igraph
graph <- graph_from_data_frame(edges, directed = FALSE)
plot(graph, vertex.size = 5, vertex.label.cex = 0.7, main = "Co-abundance Network")
