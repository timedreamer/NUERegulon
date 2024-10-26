# This script plots the log2FC comparison between maize and Arabidopsis NxTime genes.


# Output: 
# zma_version_two_species_log2FC_scatterplots_combined.pdf

# Author: Ji Huang

# 0. Prep -----------------------------------------------------------------

library(tidyverse)
library(here)
library(readxl)
library(broom)
library(cowplot)

perm_time = 1000

# 1. Load data ------------------------------------------------------------

## 1.1 Load maize and Ath shoot log2FC (+N/-N)
zma_shoot_all <- read_tsv(file = here("result", "DEgenes", "N_response", 
                                      "zma_shoot_N_repsonse_DEgenes_all.tsv.gz"))

ath_shoot_all <- read_tsv(file = here("result", "misc", "ath_result", 
                                 "ath_shoot_N_repsonse_DEgenes_all.tsv.gz"))

## 1.2 Load maize and Ath root log2FC (+N/-N)
zma_root_all <- read_tsv(file = here("result", "DEgenes", "N_response", 
                                      "zma_root_N_repsonse_DEgenes_all.tsv.gz"))

ath_root_all <- read_tsv(file = here("result", "misc", "ath_result", 
                                      "ath_root_N_repsonse_DEgenes_all.tsv.gz"))

## 1.3 Load maize and ath shoot NxTime genes. 
source(here("src", "functions", "load_all_gene-list.R"))

ath_shoot <- read_xlsx(path = here("data", "ath_time-course", 
                                   "pnas.1721487115.sd01.xlsx"), 
                       sheet = "Table S1", skip = 1)
colnames(ath_shoot) <- c("gene", "symbol", "full_name", "fdr", "JIT_cat")

## 1.4 Load maize and ath root NxTime genes. 
ath_root <- read_xlsx(path = here("data", "ath_time-course", 
                                   "pnas.1721487115.sd01.xlsx"), 
                       sheet = "Table S2", skip = 1)
colnames(ath_root) <- c("gene", "symbol", "full_name", "fdr", "JIT_cat")

## 1.5 Load ortholog table. OrthoFinder2 version.
ortho <- read_tsv(here("result", "zma_to_ath_orthogroup_long.tsv.gz")) %>% 
    dplyr::select(group, ath_gene, zma_gene) %>% distinct()

## 1.6 Load zma and ath ortholog(OF2) genes that are in the pruned network
## Shoot
output_dir <- here("result", "network_comparison", "orthofinder2_group", "shoot")
conservTF_of2 <- read_tsv(here(output_dir, "conservTF_of2.tsv"))
conservTarget_of2 <- read_tsv(here(output_dir, "conservTarget_of2.tsv"))

conserv_shoot <- unique(c(conservTF_of2$ath_gene, 
                          conservTarget_of2$ath_gene))

## Root
output_dir <- here("result", "network_comparison", "orthofinder2_group", "root")
conservTF_of2 <- read_tsv(here(output_dir, "conservTF_of2.tsv"))
conservTarget_of2 <- read_tsv(here(output_dir, "conservTarget_of2.tsv"))

conserv_root <- unique(c(conservTF_of2$ath_gene, 
                          conservTarget_of2$ath_gene))

## Load both species NxTime genes. From `008_ortholog_log2FC.R`.
load(here("result", "DEgenes", "N_response", "nxt_gene_list.RData"))

both_species_nxt_shoot |> distinct(zma_gene) |> dim()
both_species_nxt_root |> distinct(zma_gene) |> dim()

# 2. Shoot ------------------------------------------------------

pair_key_shoot <- both_species_nxt_shoot |>
    mutate(pair_key = paste0(ath_gene, "_", zma_gene)) |>
    dplyr::pull(pair_key)

## Combine ath and zma NxTime genes with log2FC.
ath_nxt_logfc_shoot <- ath_shoot_all |>
    left_join(ortho, by = c("geneid" = "ath_gene")) %>%
    drop_na(log2FoldChange, zma_gene) %>%
    dplyr::select(ath_gene = geneid,
        ath_logFC = log2FoldChange, ath_padj = padj, zma_gene) %>%
    left_join(zma_shoot_all, by = c("zma_gene" = "geneid")) %>%
    dplyr::select(starts_with("ath"), zma_gene,
        zma_logFC = log2FoldChange, zma_padj = padj) |>
    drop_na(zma_logFC)

ath_zma_log2FC <- ath_nxt_logfc_shoot |>
    filter(zma_gene %in% zma_shoot$gene) |>
    mutate(pair_key = paste0(ath_gene, "_", zma_gene)) |> 
    mutate(Type = if_else(pair_key %in% pair_key_shoot,
                          "Conserved_NxTime", "Zma_NxTime_Only"))

ath_zma_log2FC |> distinct(zma_gene, .keep_all = TRUE) |> count(Type)

plog2_shoot <- ggplot(data = ath_zma_log2FC,
    aes(x = ath_logFC, y = zma_logFC,
        color = Type, alpha = Type)) +
    geom_point() +
    ggtitle("Shoot NxTime genes") +
    scale_color_manual(values = c("green3", "grey60")) +
    scale_alpha_manual(values = c(0.4, 0.1)) +
    stat_smooth(method = "lm", se = FALSE, fullrange=TRUE) +
    xlim(c(-2.5, 7.5)) +
    ylim(c(-2.5, 7.5)) +
    coord_fixed() +
    labs(x = "Ath log2FC(+N/-N)", y = "Zma log2FC(+N/-N)") +
    cowplot::theme_cowplot() +
    theme(plot.title = element_text(size = 12))

correlations <- ath_zma_log2FC %>%
    group_by(Type) %>%
    summarize(
        correlation = cor.test(ath_logFC, zma_logFC, method = "p")$estimate,
        p_value = cor.test(ath_logFC, zma_logFC, method = "p")$p.value
    )
correlations # 0.00765 (0.681), 0.323 (p=4.26E-12)

lm(zma_logFC ~ ath_logFC, data = ath_zma_log2FC |>
    filter(Type == "Conserved_NxTime")) |> summary()

## 2. To test, whether the `both_species_nxt` are significantly better than 
## all NxT orthologs, Do a permutation test. 
all_cor <- c()
for (i in 1:perm_time) {
    set.seed(i + 585490)

    random_rep <- ath_zma_log2FC %>%
        filter(Type == "Zma_NxTime_Only") |>
        sample_n(size = length(unique(both_species_nxt_shoot$zma_gene)),
            replace = FALSE)

    logfc_cor <- tidy(cor.test(x = random_rep$ath_logFC,
        y = random_rep$zma_logFC, method = "p"))[[1, 1]]

    all_cor <- c(all_cor, logfc_cor)
}

sum(0.323<all_cor)/perm_time # p<0.001

# 3. Root ------------------------------------------------------

pair_key_root <- both_species_nxt_root |>
    mutate(pair_key = paste0(ath_gene, "_", zma_gene)) |>
    dplyr::pull(pair_key)

## Combine ath and zma NxTime genes with log2FC.
ath_nxt_logfc_root <- ath_root_all |>
    left_join(ortho, by = c("geneid" = "ath_gene")) %>%
    drop_na(log2FoldChange, zma_gene) %>%
    dplyr::select(ath_gene = geneid,
        ath_logFC = log2FoldChange, ath_padj = padj, zma_gene) %>%
    left_join(zma_root_all, by = c("zma_gene" = "geneid")) %>%
    dplyr::select(starts_with("ath"), zma_gene,
        zma_logFC = log2FoldChange, zma_padj = padj) |>
    drop_na(zma_logFC)

ath_zma_log2FC <- ath_nxt_logfc_root |>
    filter(zma_gene %in% zma_root$gene) |>
    mutate(pair_key = paste0(ath_gene, "_", zma_gene)) |> 
    mutate(Type = if_else(pair_key %in% pair_key_root,
                          "Conserved_NxTime", "Zma_NxTime_Only"))

ath_zma_log2FC |> distinct(zma_gene, .keep_all = TRUE) |> count(Type)

plog2_root <- ggplot(data = ath_zma_log2FC,
    aes(x = ath_logFC, y = zma_logFC,
        color = Type, alpha = Type)) +
    geom_point() +
    ggtitle("Root NxTime genes") +
    scale_color_manual(values = c("tan3", "grey60")) +
    scale_alpha_manual(values = c(0.4, 0.1)) +
    stat_smooth(method = "lm", se = FALSE, fullrange=TRUE) +
    xlim(c(-2, 5)) +
    ylim(c(-2, 5)) +
    coord_fixed() +
    labs(x = "Ath log2FC(+N/-N)", y = "Zma log2FC(+N/-N)") +
    cowplot::theme_cowplot() +
    theme(plot.title = element_text(size = 12))

correlations <- ath_zma_log2FC %>%
    group_by(Type) %>%
    summarize(
        correlation = cor.test(ath_logFC, zma_logFC, method = "p")$estimate,
        p_value = cor.test(ath_logFC, zma_logFC, method = "p")$p.value
    )
correlations # -0.0351 (p=0.0951), 0.471 (p=1.02E-27)

# lm(zma_logFC ~ ath_logFC, data = ath_zma_log2FC |>
#     filter(Type == "Conserved_NxTime")) |> summary()

## 2. To test, whether the `both_species_nxt` are significantly better than 
## all NxT orthologs, Do a permutation test. 
all_cor <- c()
for (i in 1:perm_time) {
    set.seed(i + 585490)

    random_rep <- ath_zma_log2FC %>%
        filter(Type == "Zma_NxTime_Only") |>
        sample_n(size = length(unique(both_species_nxt_root$zma_gene)),
            replace = FALSE)

    logfc_cor <- tidy(cor.test(x = random_rep$ath_logFC,
        y = random_rep$zma_logFC, method = "p"))[[1, 1]]

    all_cor <- c(all_cor, logfc_cor)
}

sum(0.471 < all_cor)/perm_time # p<0.001


# 4. Save plots -----------------------------------------------------------

pdf(here("result", "DEgenes", "N_response",
         "zma_version_two_species_log2FC_scatterplots_combined.pdf"), 
    width = 5, height = 3.5)

plog2_shoot
plog2_root
dev.off()