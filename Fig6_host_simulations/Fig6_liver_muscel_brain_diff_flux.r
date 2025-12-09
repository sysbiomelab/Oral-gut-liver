#########################################################
#################### liver model ########################
#########################################################

library(readxl)
cirrhosis_Ex <- read_excel("~/.../Ex.xlsx", 
                           sheet = "LC")

normal_Ex <- read_excel("~/.../Ex.xlsx", 
                        sheet = "normal")

library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)


eps <- 1e-6 
thr <- 3.5 
set.seed(1)

num_cols <- c("0","5")  
clean_numeric <- function(df) {
  df %>%
    mutate(across(all_of(num_cols),
                  ~ suppressWarnings(as.numeric(str_replace_all(., ",", "")))))
}

extract_metabolite <- function(eqs) {
  eqs %>%
    str_remove("<=>.*$") %>%
    str_remove("\\[.*$") %>%
    str_squish()
}

normal  <- clean_numeric(normal_Ex)  %>% mutate(metabolite = extract_metabolite(eqs))
cirrh   <- clean_numeric(cirrhosis_Ex) %>% mutate(metabolite = extract_metabolite(eqs))

#################### log2FC (cirrhosis vs normal) ########################
merged_cn <- normal %>%
  select(rxns, metabolite, n0 = `0`, n5 = `5`) %>%
  inner_join(
    cirrh %>% select(rxns, c0 = `0`, c5 = `5`),
    by = "rxns"
  ) %>%
  mutate(
    log2FC_0 = log2(c0 + eps) - log2(n0 + eps),   # ammonia = 0
    log2FC_5 = log2(c5 + eps) - log2(n5 + eps)    # ammonia = 5
  )

merged_cn <- merged_cn %>%
  mutate(label_flag = abs(log2FC_0) >= thr | abs(log2FC_5) >= thr)

lab_df <- merged_cn %>% filter(label_flag)

# x: ammonia=5 log2FC，y: ammonia=0 log2FC
p <- ggplot(merged_cn, aes(x = log2FC_5, y = log2FC_0)) +
  geom_point(aes(color = label_flag),
             size = 2.5, alpha = 0.9,
             position = position_jitter(width = 0.05, height = 0.05)) +
  ggrepel::geom_text_repel(
    data = lab_df,
    aes(label = metabolite),
    size = 2, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  geom_hline(yintercept = 0, linetype = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  coord_equal() +
  scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "tomato"), guide = "none") +
  labs(
    x = "log2FC at ammonia = 5  (cirrhosis vs normal)",
    y = "log2FC at ammonia = 0  (cirrhosis vs normal)",
    subtitle = paste0("Labeled if |log2FC| ≥ ", thr)
  ) +
  theme_minimal(base_size = 12)

print(p)

write.csv(merged_cn,"liver_GSMM_2_vs_0_stat.csv")


#########################################################
################ muscle and brain model #################
#########################################################
library(readxl)
muscle_Ex <- read_excel("~/.../muscle_bran_Ex.xlsx", 
                        sheet = "muscle")
brain_Ex <- read_excel("~/.../muscle_bran_Ex.xlsx", 
                       sheet = "brain")
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

num_cols <- as.character(1:12)
clean_numeric <- function(df) {
  df %>%
    mutate(across(all_of(num_cols),
                  ~ suppressWarnings(as.numeric(str_replace_all(., ",", "")))))
}

muscle <- clean_numeric(muscle_Ex)
brain  <- clean_numeric(brain_Ex)

# ---  log2FC(12 vs 1) ---
eps <- 1e-6
muscle <- muscle %>%
  mutate(log2FC_muscle = log2(`12` + eps) - log2(`1` + eps))
brain  <- brain %>%
  mutate(log2FC_brain  = log2(`12` + eps) - log2(`1` + eps))

extract_metabolite <- function(eqs) {
  eqs %>%
    str_remove("<=>.*$") %>%
    str_remove("\\[.*$") %>%
    str_squish()
}
muscle <- muscle %>% mutate(metabolite = extract_metabolite(eqs))
brain  <- brain  %>% mutate(metabolite = extract_metabolite(eqs))

# --- join（by rxns）---
merged <- muscle %>%
  select(rxns, metabolite, log2FC_muscle) %>%
  inner_join(brain %>% select(rxns, log2FC_brain), by = "rxns") %>%
  filter(is.finite(log2FC_muscle), is.finite(log2FC_brain))


library(dplyr)
library(ggplot2)
library(ggrepel)

thr <- 1  # |log2FC| threshold

# 1) label
merged2 <- merged %>%
  mutate(label_flag = abs(log2FC_brain) >= thr | abs(log2FC_muscle) >= thr)

# 2) label 
lab_df <- merged2 %>% filter(label_flag)


############
# overlapped avoid
set.seed(1)  
ggplot(merged2, aes(x = log2FC_brain, y = log2FC_muscle)) +
  geom_point(
    aes(color = label_flag),
    size = 3, alpha = 0.9,
    position = position_jitter(width = 0.05, height = 0.05)
  ) +
  ggrepel::geom_text_repel(
    data = lab_df,
    aes(label = metabolite),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  coord_equal() +
  scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "tomato"), guide = "none") +
  labs(x = "Brain log2FC (100 vs 0)", y = "Muscle log2FC (100 vs 0)",
       # title = "Metabolite log2FC: Brain, Muscle",
       subtitle = paste0("Labeled if |log2FC| ≥ ", thr))+
  theme_minimal(base_size = 12)
#  volcano_brain_muscle 4*4

write.csv(merged2,"brain_muscle_GSMM_stat.csv")


#############
# filter the non-zero ones on both side, 100 and 0
library(dplyr)
library(ggplot2)
library(ggrepel)

thr <- 1

# 
merged_full <- muscle %>%
  select(rxns, metabolite, m1 = `1`, m12 = `12`, log2FC_muscle) %>%
  inner_join(
    brain %>% select(rxns, b1 = `1`, b12 = `12`, log2FC_brain),
    by = "rxns"
  )

merged_nz <- merged_full %>%
  filter(is.finite(m1), is.finite(m12), is.finite(b1), is.finite(b12)) %>%
  filter(m1 != 0, m12 != 0, b1 != 0, b12 != 0)

merged_nz <- merged_nz %>%
  mutate(label_flag = abs(log2FC_brain) >= thr | abs(log2FC_muscle) >= thr)

lab_df <- merged_nz %>% filter(label_flag)

set.seed(1)  
ggplot(merged_nz, aes(x = log2FC_brain, y = log2FC_muscle)) +
  geom_point(
    aes(color = label_flag),
    size = 2.2, alpha = 0.9,
    position = position_jitter(width = 0.05, height = 0.05)
  ) +
  ggrepel::geom_text_repel(
    data = lab_df,
    aes(label = metabolite),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  coord_equal() +
  scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "tomato"), guide = "none") +
  labs(
    x = "Brain log2FC (100 vs 0)",
    y = "Muscle log2FC (100 vs 0)",
    subtitle = paste0("labeled if |log2FC| ≥ ", thr)
  ) +
  theme_minimal(base_size = 12)

