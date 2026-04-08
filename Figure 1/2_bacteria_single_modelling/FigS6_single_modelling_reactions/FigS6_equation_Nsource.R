library(readxl)
NH3_pathway_table <- read_excel("~/.../NH3_pathway_table.xlsx", 
                                sheet = "data")


# Load necessary libraries
library(tidyverse)
library(ggplot2)


# Ensure numeric columns are selected for reshaping
NH3_pathway_long <- NH3_pathway_table %>%
  select(Equation, N_source, where(is.numeric)) %>%
  pivot_longer(cols = -c(Equation, N_source), # Exclude 'Equation' and 'N_source'
               names_to = "Microbe",
               values_to = "Value")

# View the reshaped data
head(NH3_pathway_long)

NH3_pathway_long2 <- NH3_pathway_long %>%
  mutate(Value = case_when(
    is.na(Value) ~ 0,                      # If NA, set to 0
    Value >= 0.1 ~ 0.1,                    # If value >= 0.1, set to 0.1
    Value >= 0.01 & Value < 0.1 ~ 0.01,    # If 0.01 <= value < 0.1, set to 0.01
    Value >= 0.001 & Value < 0.01 ~ 0.001, # If 0.001 <= value < 0.01, set to 0.001
    Value > 0 & Value < 0.001 ~ 0.0001,    # If 0 < value < 0.001, set to 0.0001
    TRUE ~ Value                           # If none of the conditions are met, keep original value
  ))


# Step 2: Visualizing the Data
# Heatmap
ggplot(NH3_pathway_long2, aes(x = Equation, y = Microbe, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Heatmap of NH3 Pathways Across Microbial Species",
       x = "Pathway",
       y = "Microbial Species",
       fill = "Nitrogen Value")

############


NH3_pathway_long3 <- NH3_pathway_long %>%
  mutate(Value_category = case_when(
    is.na(Value) ~ "None",  # If the value is NA, assign "None"
    TRUE ~ cut(Value,
               breaks = c(0, 0.001, 0.01, 0.1, 0.5, Inf), 
               labels = c("Very Low", "Low", "Medium", "High", "Very High"),
               include.lowest = TRUE)
  ))
# Now you can use `Value_category` for categorical visualization
colorRampPalette(c("#A6D0FF", "#005B96"))(5)  # Light blue to dark blue
c("#A6D0FF", "#7CB2E4", "#5395CA", "#2978B0", "#005B96")
ggplot(NH3_pathway_long3, aes(x = Equation, y = Microbe, fill = Value_category)) +
  geom_tile() +
  # scale_fill_manual(values = c("Very Low" = "red", "Low" = "orange", 
  #                              "Medium" = "yellow", "High" = "green", "Very High" = "blue","None" = "grey")) +
  scale_fill_manual(values = c("Very High" = "#005B96",  # Dark Forest Green
                               "High" = "#2978B0",     # Dark Green
                               "Medium" = "#5395CA",   # Medium Green
                               "Low" = "#7CB2E4",      # Light Olive Green
                               "Very Low" = "#A6D0FF",  # Light Green
                               "None" = "grey"), 
                    limits = c("Very High", "High", "Medium", "Low", "Very Low", "None")) +  # Custom legend order
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "NH3 Pathways Across Microbial Species",
       x = "Pathway", 
       y = "Microbial Species", 
       fill = "Category")+
  coord_flip()
# save pathway_ggheatmap 10*7

# Install the necessary package if not already installed
# install.packages("pheatmap")

# Load the necessary libraries
# Load necessary libraries
library(tidyverse)
library(pheatmap)


# Step 1: Categorize the 'Value' column using 'cut()' into specified categories
NH3_pathway_long3

# Step 2: Reshape the data to have bacteria as rows and pathways as columns
# Here we pivot the data so that rows are bacteria, columns are pathways, and values are the categories

NH3_pathway_wide <- NH3_pathway_long3 %>%
  select(Microbe, Equation, Value_category) %>% # Keep Microbe, Pathway, and Value_category
  pivot_wider(names_from = Equation, values_from = Value_category, values_fn = list(Value_category = first)) 

# Step 3: Convert the data into a matrix for the heatmap (we need it in numeric format for pheatmap)
# Convert categorical values to numeric values: "Very Low" -> 1, "Low" -> 2, etc.
value_mapping <- c("None" = 0,"Very Low" = 1, "Low" = 2, "Medium" = 3, "High" = 4, "Very High" = 5)

NH3_pathway_matrix <- NH3_pathway_wide %>%
  select(-Microbe) %>%
  mutate(across(everything(), ~ value_mapping[.])) %>%
  as.matrix()


# Set row names of the matrix to be the 'Microbe' column values
rownames(NH3_pathway_matrix) <- NH3_pathway_wide$Microbe

######-------------- heatmap --------------######
# Define the legend breaks and labels based on your value_mapping
legend_breaks <- c(0, 1, 2, 3, 4, 5)
legend_labels <- c("None", "Very Low", "Low", "Medium", "High", "Very High")

# Create a blue-based color palette
blue_palette <- colorRampPalette(c("#A6D0FF", "#005B96"))(5)  # Light blue to dark blue
c("#A6D0FF", "#7CB2E4", "#5395CA", "#2978B0", "#005B96")

# Step 4: Generate the clustered heatmap using pheatmap
pheatmap(t(NH3_pathway_matrix), 
         clustering_distance_rows = "euclidean",  # Row clustering method (Euclidean distance)
         clustering_distance_cols = "euclidean",  # Column clustering method (Euclidean distance)
         clustering_method = "complete",           # Clustering method (complete linkage)
         scale = "none",                          # No scaling, but you can scale if needed
         color = colorRampPalette(c("lightgrey","#A6D0FF", "#7CB2E4", "#5395CA", "#2978B0", "#005B96"))(6),  # Green color gradient
         # main = "NH3 Pathways Across Microbial Species",
         show_rownames = TRUE,                    # Show row names (bacteria)
         show_colnames = TRUE,                    # Show column names (pathways)
         legend = TRUE,                           # Show color legend
         # breaks = seq(0, 5, by = 1),
         legend_breaks = legend_breaks,
         legend_labels = legend_labels, # Define color breaks for categorical data
         treeheight_row = 10,  # Disable row dendrogram
         treeheight_col = 10,  # Disable column dendrogram
         cellwidth = 20,  # Adjust width of the cells
         cellheight = 20) # Adjust height of the cells              
# save Pathway_pheatmap 11*7.2 L

###### add original N source ####

