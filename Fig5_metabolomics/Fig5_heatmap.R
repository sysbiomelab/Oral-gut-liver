library(readxl)
sample_info <- read_excel("~/.../stool_samples.xlsx", 
                          sheet = "Sheet1")

library(readxl)
blood_metabolomics <- read_excel("~/.../Metabolomics Portal vein TIPS.xlsx")

library(readxl)
metaG_16 <- read_excel("~/.../Kraken2_Data.xlsx", 
                       sheet = "Sheet1")
metaG_16 = metaG_16[-1,]

library(dplyr)
library(stringr)
library(tidyr)

## --- 0) clean types ---------------------------------------------------
metaG_clean <- metaG_16 %>%
  mutate(across(-Sample, ~ suppressWarnings(as.numeric(as.character(.)))))

blood_clean <- blood_metabolomics %>%
  mutate(Patient = as.character(Patient)) %>%
  mutate(across(-Patient, ~ suppressWarnings(as.numeric(as.character(.)))))

## --- 1) relative abundance for metaG ----------------------------------
# locate the total row
bact_row <- which(metaG_clean$Sample == "Bacteria <bacteria>")
stopifnot(length(bact_row) == 1)

# per-column totals
totals <- as.numeric(metaG_clean[bact_row, -1])

# avoid division-by-zero
totals[totals == 0] <- NA

# remove the total row and divide every column by its total
metaG_rel <- metaG_clean[-bact_row, ]
metaG_rel[, -1] <- sweep(metaG_rel[, -1, drop = FALSE], 2, totals, "/")
# (optional) to get percentages: metaG_rel[, -1] <- metaG_rel[, -1] * 100

## --- 2) build mapping: column -> patient & replicate ------------------
library(dplyr)
library(stringr)

metaG_cols <- setdiff(colnames(metaG_rel), "Sample")

# 
m <- str_match(metaG_cols, "PRE_(\\d+)_(\\d+)")   # col × (full, grp1, grp2)

col_map <- tibble(
  metaG_col = metaG_cols,
  MGP_ID = paste0("PRE_", m[,2]),                 
  rep    = as.integer(m[,3])                      
) %>%
  left_join(sample_info %>% rename(Sample_ID = `Sample ID`),
            by = c("MGP_ID" = "MGP ID")) %>%
  mutate(obs_id = paste0("P", Sample_ID, "_r", rep))


## --- 3) metaG matrix: features × observations (patient-replicates) ----
metaG_obs <- metaG_rel %>%
  select(Sample, all_of(metaG_cols))
colnames(metaG_obs)[-1] <- col_map$obs_id

G_mat <- as.matrix(metaG_obs[, -1])
rownames(G_mat) <- metaG_obs$Sample

## --- 4) align patients ------------------------------------------


# patient with paired blood and stool samples
patients_blood  <- setdiff(names(blood_wide), "metabolite")
patients_metaG  <- as.character(unique(col_map$Sample_ID))
patients_common <- intersect(patients_blood, patients_metaG)

# missing paired samples patients
missing_in_blood <- setdiff(patients_metaG, patients_blood)
missing_in_metaG <- setdiff(patients_blood, patients_metaG)
if(length(missing_in_blood)) message("No blood for patients: ", paste(missing_in_blood, collapse=", "))
if(length(missing_in_metaG)) message("No metaG for patients: ", paste(missing_in_metaG, collapse=", "))

# 
col_map_keep <- col_map %>%
  filter(as.character(Sample_ID) %in% patients_common)

# 
metaG_cols_keep <- col_map_keep$metaG_col
metaG_obs <- metaG_rel %>%
  select(Sample, all_of(metaG_cols_keep))
colnames(metaG_obs)[-1] <- col_map_keep$obs_id

G_mat <- as.matrix(metaG_obs[, -1, drop = FALSE])
rownames(G_mat) <- metaG_obs$Sample

# order: col_map_keep 
obs_order <- col_map_keep$obs_id
obs_pat   <- as.character(col_map_keep$Sample_ID)

B_mat <- blood_wide %>%
  select(metabolite, all_of(patients_common)) %>%
  {
    mm <- sapply(seq_along(obs_order), function(k){
      patient_k <- obs_pat[k]
      as.numeric(.[[patient_k]])
    })
    mm <- matrix(mm, nrow = nrow(.), byrow = FALSE)
    rownames(mm) <- .$metabolite
    colnames(mm) <- obs_order
    mm
  }

## --- 5) Spearman  -----------------------------------------------

common_obs <- intersect(colnames(B_mat), colnames(G_mat))
B <- B_mat[, common_obs, drop = FALSE]
G <- G_mat[, common_obs, drop = FALSE]

spearman_pairwise <- function(X, Y){
  rho  <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(Y))
  pval <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(Y))
  for(i in seq_len(nrow(X))){
    x <- as.numeric(X[i, ])
    for(j in seq_len(nrow(Y))){
      y <- as.numeric(Y[j, ])
      ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      rho[i, j]  <- unname(ct$estimate)
      pval[i, j] <- ct$p.value
    }
  }
  rownames(rho)  <- rownames(X)  # metabolites
  colnames(rho)  <- rownames(Y)  # metaG features
  rownames(pval) <- rownames(X)
  colnames(pval) <- rownames(Y)
  list(rho = rho, p = pval)
}

res <- spearman_pairwise(B, G)
rho_mat <- res$rho
p_mat   <- res$p

# BH adj
p_adj_mat <- t(apply(p_mat, 1, function(v) p.adjust(v, method = "BH")))

# write.csv(rho_mat,   "rho_spearman_relab_replicates.csv")
# write.csv(p_mat,     "pval_spearman_relab_replicates.csv")
# write.csv(p_adj_mat, "pval_BH_spearman_relab_replicates.csv")


###################-------------------------------------######################################
# filter
keep_idx <- !grepl("^M\\d+T\\d+", rownames(rho_mat))

rho_mat_annot <- rho_mat[keep_idx, , drop = FALSE]
p_mat_annot   <- p_mat[keep_idx, , drop = FALSE]

# FDR adj
p_adj_mat_annot <- t(apply(p_mat_annot, 1, function(v) p.adjust(v, method = "BH")))

# save
# write.csv(rho_mat_annot,   "rho_spearman_annotated.csv")
# write.csv(p_mat_annot,     "pval_spearman_annotated.csv")
# write.csv(p_adj_mat_annot, "pval_BH_spearman_annotated.csv")



### plot ###
### plot ###
### plot ###
### plot ###
### plot ###
### plot ###


##########--------------------------------------------------------------#########
### plot ###
library(pheatmap)

## ---- 1) keep annotated metabolites ------------------------------------------
# # rm MxxxTxxx rows
# keep_idx <- !grepl("^M\\d+T\\d+", rownames(rho_mat))
# rho_mat_annot <- rho_mat[keep_idx, , drop = FALSE]
# p_mat_annot   <- p_mat[keep_idx, , drop = FALSE]

## ---- 2) keep p < 0.05 metabolites -------------------------------------------
sig_idx <- apply(p_mat_annot, 1, function(x) any(x < 0.05, na.rm = TRUE))
rho_sig <- rho_mat_annot[sig_idx, , drop = FALSE]
p_sig   <- p_mat_annot[sig_idx, , drop = FALSE]



## ---- 3) generate significance matrix ----------------------------------------
signif_mat <- matrix("", nrow = nrow(p_sig), ncol = ncol(p_sig))
signif_mat[p_sig < 0.05]  <- "*"
signif_mat[p_sig < 0.01]  <- "**"
signif_mat[p_sig < 0.001] <- "***"
rownames(signif_mat) <- rownames(p_sig)
colnames(signif_mat) <- colnames(p_sig)


## ---- 4) plot heatmap --------------------------------------------

my_colors <- colorRampPalette(c("darkblue", "white", "darkred"))(200)

# set breaks
my_breaks <- seq(-0.5, 0.5, length.out = 201)

pheatmap(
  rho_sig,
  cluster_rows = TRUE,    
  cluster_cols = TRUE,    
  display_numbers = signif_mat, 
  number_color = "black", 
  color = my_colors,
  breaks = my_breaks,
  main = "Spearman correlation (annotated metabolites, p<0.05)",
  cellwidth  = 12,   
  cellheight = 12    
)



###################-------------plotting-------------------################
library(readxl)
prevalence_translocation_species_TIPS <- read_excel("~/.../prevalence_translocation_species_TIPS.xlsx")


library(ComplexHeatmap)
library(circlize)
library(grid)

library(ComplexHeatmap)
library(circlize)
library(grid)

## ---- 1) add prevalence data --------------------------------------
# keep order consistant with rho_sig column name
species_order <- colnames(rho_sig)
prev_vec <- prevalence_translocation_species_TIPS$prevalence_overall[
  match(species_order, prevalence_translocation_species_TIPS$Group)
]
names(prev_vec) <- species_order

## ---- 2) colour（0=white） ----------------------------------
min_val <- min(rho_sig, na.rm = TRUE)
max_val <- max(rho_sig, na.rm = TRUE)

col_fun <- colorRamp2(c(min_val, 0, max_val), c("blue", "white", "red"))

## ---- 3) label significance ------------------------------------------
signif_mat <- matrix("", nrow = nrow(p_sig), ncol = ncol(p_sig))
signif_mat[p_sig < 0.05]  <- "*"
signif_mat[p_sig < 0.01]  <- "**"
signif_mat[p_sig < 0.001] <- "***"
rownames(signif_mat) <- rownames(p_sig)
colnames(signif_mat) <- colnames(p_sig)

## ---- 4) anno：barplot -------------------------------------

ha_col <- HeatmapAnnotation(
  "Prevalence in ACLD" = anno_barplot(
    prev_vec,
    border = FALSE,
    gp = gpar(fill = "grey40"),
    height = unit(2, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 12)  # smaller font
)


## ---- 5)  heatmap ----------------------------------------------

ht <- Heatmap(
  rho_sig,
  name = "rho",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = ha_col,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(signif_mat[i, j], x, y, gp = gpar(fontsize = 7))
  },
  width = unit(ncol(rho_sig)*5, "mm"),   # ensure width proportional to ncol
  height = unit(nrow(rho_sig)*5, "mm")   # ensure height proportional to nrow
)

draw(ht)

#######------------------------------######








