# This script plotted the PCC increase using random 23/24 genes versus
# MYB-DIV conserved genes.
# It also did the min-max normalization for the XGboost results.

# I plotted the maize and Arabidopsis version. I chose the Okabe&Ito color palette.
# Maize: #E69F00; Ath: #009773

# Author: Ji Huang
# Date: 2024-03-07

# Output:
# dumbbell_2species_MYBDIVgenes_with5Random.pdf
# MYBDIV_model_norm_XGvalues.tsv

# 0. Prep -----------------------------------------------------------------

library(here)
library(tidyverse)
library(ggalt)
library(cowplot)
library(purrr)

outdir = here("result", "machine_learning", "xgboost")

new_max = 10
new_min = 1

# 1. Maize dumbbell plot ---------------------------------------------------

## 1.1 Load conserved genes result
load(file = here(outdir,
                 "XGBoost.Zm-leaf-orthoMYBDIV-TotalNUE-output.RData"))

genotype <- unlist(strsplit(rownames(data)[seq(2, by = 2, 32)], "_", fixed = TRUE
))[seq(2, by = 3, 48)]

genotype = factor(genotype, levels = c("B73xIHP1", "B73xILP1", "B73xLH82", "B73xMo17",
                                       "B73xMo18W", "B73xOh7B", "B73xPH207", "B73xPHG47", "B73xPHG84",
                                       "B73", "IHP1", "ILP1", "LH82", "Mo17", "PH207", "PHG84"))

df24 <- data.frame(genotype = genotype,
                   Correlation = COR)

## 1.2 Load five random genes result
random_seeds <- c(165, 524, 1532, 1746, 2024)

random_cors <- map(random_seeds, function(seed) {
  load(file = here(outdir, paste0("XGBoost.Zm-leaf-MYBDIVrandom24_", seed, "-TotalNUE-output.RData")))
  saved_items$COR
})

names(random_cors) <- paste0("random", 1:5)

df24 <- df24 %>%
  bind_cols(as.data.frame(random_cors)) %>%
  mutate(Random24Genes = rowMeans(select(., starts_with("random")))) %>%
  dplyr::rename(ConservedModule = Correlation)

## 1.3 Plot
p_zma <- ggplot(df24, aes(y = genotype, x = Random24Genes, xend = ConservedModule)) +
  ggalt::geom_dumbbell(size = 1.5, color = "#e3e2e1",
                       size_x = 3.2, size_xend = 3.2,
                       colour_x = "#ffc034", colour_xend = "#cd8d00",
                       dot_guide = TRUE, dot_guide_size = 0.1) +
  geom_point(aes(x = random1), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random2), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random3), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random4), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random5), color = "grey60", alpha = 0.4, size = 1.5) +
  labs(x = "Correlation Coefficient", y = "Zma Genotype") +
  theme_cowplot()

## Wilcox test for the p-value
wilcox.test(df24$ConservedModule, df24$Random24Genes, paired = TRUE, alternative= c("greater")) 
# V = 136, p-value = 1.526e-05


# 2. Arabidopsis dumbbell plot ---------------------------------------------

## 2.1 Load conserved genes result
load(file = here(outdir,
    "XGBoost.Ath-leaf-orthoMYBDIV-NUE-output.RData"))

genotype <- factor(unlist(strsplit(rownames(data)[1:18], "_", fixed = T
))[seq(1, by = 3, 54)])

df23 <- data.frame(genotype = genotype,
    Correlation = COR)

## 2.2 Load five random genes result 
random_seeds <- c(1404, 1931, 525, 913, 117)

random_cors <- map(random_seeds, function(seed) {
  load(file = here(outdir, paste0("XGBoost.Ath-leaf-random23_bothNxT_", seed, "-TotalNUE-output.RData")))
  saved_items$COR
})

names(random_cors) <- paste0("random", 1:5)

df23 <- df23 %>%
  bind_cols(as.data.frame(random_cors)) %>%
  mutate(Random23Genes = rowMeans(select(., starts_with("random")))) %>%
  dplyr::rename(ConservedModule = Correlation)

## 2.3 Plot
p_ath <- ggplot(df23, aes(y = genotype, x = Random23Genes, xend = ConservedModule)) +
  ggalt::geom_dumbbell(size = 1.5, color = "#e3e2e1",
                       size_x = 3.2, size_xend = 3.2,
                       colour_x = "#00d198", colour_xend = "#008560",
                       dot_guide = TRUE, dot_guide_size = 0.1) +
  geom_point(aes(x = random1), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random2), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random3), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random4), color = "grey60", alpha = 0.4, size = 1.5) +
  geom_point(aes(x = random5), color = "grey60", alpha = 0.4, size = 1.5) +
  labs(x = "Correlation Coefficient", y = "Ath Genotype") +
  theme_cowplot()

## Wilcox test for the p-value
wilcox.test(df23$ConservedModule, df23$Random23Genes, paired = TRUE, alternative = c("greater")) 
# V = 154, p-value = 0.0007896

# 3. Min-Max normalize the Zma importance score -------------------------------

zma_score <- read_tsv(here("result", "machine_learning", "xgboost",
    "orthoMYBDIV_ZmTotalNUE.XGBoost-importantgene-Athomolog-frequency.tsv"),
col_names = c("Gene", "XGvalue", "gene_name", "gene_anno", "ath_ortho", "N")
)

zma_score <- zma_score %>%
    select(-ath_ortho, -N) %>%
    mutate(norm_XGvalue = (XGvalue - min(XGvalue)) / (max(XGvalue) - min(XGvalue)) * (new_max - new_min) + new_min)

# 4. Min-Max normalize the Arabidopsis importance score -------------------------------

ath_score <- read_tsv(here("result", "machine_learning", "xgboost",
    "orthoMYBDIV_23gene_AthNUE.XGBoost-importantgene-Athomolog-frequency.tsv"),
col_names = c("Gene", "XGvalue", "gene_name", "gene_anno", "N"),
skip = 1
)

ath_score <- ath_score %>%
    select(-N) %>%
    mutate(norm_XGvalue = (XGvalue - min(XGvalue)) / (max(XGvalue) - min(XGvalue)) * (new_max - new_min) + new_min)


# 5. Combine and save -----------------------------------------------------

# 5.1 Combine the plots
combined_plot <- plot_grid(p_zma, p_ath, align = "h")

ggsave(here(outdir, "dumbbell_2species_MYBDIVgenes_with5Random.pdf"),
    combined_plot, width = 8, height = 4)


# 5.2 Combine XGboost values
combine_score <- bind_rows(zma_score, ath_score)

write_tsv(combine_score, here(outdir, "MYBDIV_model_norm_XGvalues.tsv"))
