############################################################
## Dirichlet-multinomial reactotype analysis
## (for gut or oral reactobiome profiles)
############################################################

## 1. Load reactobiome count table -------------------------
## The input table should have:
##   - rows = features (e.g. reactions)
##   - columns = samples
## Replace the path with your own file.
library(readr)

reactobiome_mat <- read.csv(
  "xxx/path_to_data/reactobiome_frommsp_gut_or_oral.csv",
  row.names = 1,
  check.names = FALSE
)

## DirichletMultinomial expects samples in rows and features in columns,
## so we transpose the matrix:
count_matrix <- t(reactobiome_mat)

############################################################
## 2. Fit DMN models with different k and select best model
############################################################

library(DirichletMultinomial)
library(parallel)

# Fit DMN models for k = 1 to 10 components
set.seed(123)
Fit_list <- mclapply(
  1:10,
  dmn,
  count   = count_matrix,
  verbose = TRUE
)

# Compute Laplace approximation of the log-likelihood for model comparison
lplc <- sapply(Fit_list, laplace)

# Select the model with the best (lowest) Laplace value
best_index <- which.min(lplc)
best_model <- Fit_list[[best_index]]

cat("Selected number of reactotypes (k):", best_index, "\n")

# Plot Laplace values across components
plot(
  lplc,
  type = "b",
  xlab = "Number of Dirichlet Components (k)",
  ylab = "Model Fit (Laplace)"
)

############################################################
## 3. Compare fitted probabilities (k = 1 vs best model)
############################################################

# Fitted component probabilities for k = 1 (null model) and best model
p0 <- fitted(Fit_list[[1]],      scale = TRUE)   # 1-component model
pK <- fitted(best_model,         scale = TRUE)   # best model
colnames(pK) <- paste0("m", 1:ncol(pK))

# Difference between 1-component model and best model
meandiff <- colSums(abs(pK - as.vector(p0)))     # per component
diff     <- rowSums(abs(pK - as.vector(p0)))     # per feature
o        <- order(diff, decreasing = TRUE)
cdiff    <- cumsum(diff[o]) / sum(diff)

# Top 40 features contributing most to between-model differences
df_top40 <- head(
  cbind(
    Mean  = p0[o],
    pK[o, ],
    diff  = diff[o],
    cdiff = cdiff
  ),
  40
)

# Full table for all features (if needed)
df_full <- cbind(
  Mean  = p0[o],
  pK[o, ],
  diff  = diff[o],
  cdiff = cdiff
)

# Optional: write the tables to disk
# write.csv(df_top40, "xxx/path_to_output/reactotype_diff_top40.csv")
# write.csv(df_full,  "xxx/path_to_output/reactotype_diff_full.csv")

# Optional: if a custom heatmap function is available, it can be used as:
# heatmapdmn(count_matrix, Fit_list[[1]], best_model, n_features = 40, lblwidth = 60)

############################################################
## 4. Assign samples to reactotypes (clusters)
############################################################

# Posterior group membership probabilities for each sample
# Rows = samples, columns = components
post_prob <- best_model@group

# Assign each sample to the component (reactotype) with maximum posterior probability
reactotype_assignment <- apply(post_prob, 1, which.max)

reactotype_df <- data.frame(
  SampleID   = rownames(post_prob),
  reactotype = factor(reactotype_assignment)
)

# Optional: save reactotype labels
# write.csv(reactotype_df,
#           "xxx/path_to_output/reactotype_assignment.csv",
#           row.names = FALSE)
