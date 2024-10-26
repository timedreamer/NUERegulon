# This scirpt identifies JIT genes.

# Output:
# maize_time_course\result\JIT\JIT_genes\

# Auhtor: Ji Huang
# Date: 2022-01-28

# 0. Prep -----------------------------------------------------------------

library(tidyverse)
library(here)
library(ggbeeswarm)
theme_set(cowplot::theme_cowplot())

source(here("src", "functions", "load_read-count_design-table.R"))

FC_cutoff <- 2

output_dir <- here("result", "JIT", "JIT_genes")

# 1. Calculate e Fold Change -----------------------------------------------

## First need to load the quant normalized data and calculate log2FC
mts_quant <- read_tsv(here("result", "quant_normalized_expr_all.tsv.gz"))

nxt_gene <- read_tsv(here("result", "DEgenes", "edge", "zma_NxTgenes_df5_edge.tsv"))

mts_quant_nxt <- mts_quant %>% filter(geneid %in% nxt_gene$gene)

mts_quant_nxt <- mts_quant_nxt %>% select(sort(colnames(mts_quant_nxt))) %>% select(-(1:3))


mts_quant_nxt <- mts_quant_nxt %>% 
    pivot_longer(col=-geneid, names_to = "condition", values_to = "quant") %>% 
    separate(condition, into = c("condition", "Ntreatment"), sep="-") %>% 
    separate(condition, sep = "mS", into = c("time", "rep"))

## calculate FC.
nxt_foldChange <- mts_quant_nxt %>% group_by(geneid, time, Ntreatment) %>% 
    summarise(mean_quant = mean(quant)) %>% 
    ungroup() %>% 
    group_by(geneid, time) %>%  # these two may not be necessary.
    mutate(change=((mean_quant+1) /(mean_quant[Ntreatment=="noN"]+1))) %>% 
    ungroup() %>% 
    filter(Ntreatment == "N") %>% 
    select(-Ntreatment, -mean_quant)


nxt_foldChange <- nxt_foldChange %>% 
    mutate(log2FC = log2(change)) %>% 
    mutate(time = as.numeric(time)) %>% 
    select(-change) %>% 
    filter(abs(log2FC) > abs(log2(FC_cutoff)))

# 2. Bin genes into JIT bins ----------------------------------------------

mts_jit <- nxt_foldChange %>% group_by(time) %>% 
    group_split()

names(mts_jit) <- nxt_foldChange %>% group_by(time) %>% 
    group_keys() %>% mutate(time = paste0("g", time)) %>% pull(time)

## Bin each gene into a time slot.

jit_5 <- mts_jit$g5$geneid
used_gene <- jit_5

jit_10 <- setdiff(mts_jit$g10$geneid, used_gene) 
used_gene <- c(used_gene, jit_10)

jit_15 <- setdiff(mts_jit$g15$geneid, used_gene) 
used_gene <- c(used_gene, jit_15)

jit_20 <- setdiff(mts_jit$g20$geneid, used_gene) 
used_gene <- c(used_gene, jit_20)

jit_30 <- setdiff(mts_jit$g30$geneid, used_gene) 
used_gene <- c(used_gene, jit_30)

jit_45 <- setdiff(mts_jit$g45$geneid, used_gene) 
used_gene <- c(used_gene, jit_45)

jit_60 <- setdiff(mts_jit$g60$geneid, used_gene) 
used_gene <- c(used_gene, jit_60)

jit_90 <- setdiff(mts_jit$g90$geneid, used_gene) 
used_gene <- c(used_gene, jit_90)

jit_120 <- setdiff(mts_jit$g120$geneid, used_gene) 
used_gene <- c(used_gene, jit_120)


# 3. Save JIT result ------------------------------------------------------

# Convert the lists to tibble. Remove lncRNAs.
jit_all <- jit_all %>% enframe(name = "time", value = "gene") %>% 
    unnest(cols = "gene") %>% 
    filter(!grepl("^E", .$gene))

jit_all %>% dplyr::count(time)

## Save JIT all result.
write_tsv(jit_all, file = here(output_dir, 
                               paste0("maize_JIT_edge_all_genes_FCcutoff",
                                      FC_cutoff, ".tsv"))
          )

# 4. Save JIT Fold-change values ------------------------------------------

jit_shoot <- read_tsv(here(output_dir, "maize_JIT_edge_all_genes_FCcutoff2.tsv")) %>% 
    mutate(time = factor(time, levels = c("j5", "j10", "j15", "j20", 
                                          "j30", "j45", "j60", 
                                          "j90", "j120")))

jit_shoot <- jit_all
nxt_foldChange_1 <- nxt_foldChange %>% 
    mutate(time_cat = paste0("j", time)) %>% 
    mutate(gene_time_id = paste0(geneid, "_", time_cat)) %>% 
    select(gene_time_id, log2FC)


jit_shoot_fc <- jit_shoot %>% 
    mutate(gene_time_id = paste0(gene, "_", time)) %>% 
    left_join(nxt_foldChange_1, by = "gene_time_id") 

## Save result.
write_tsv(jit_shoot_fc, file = here(output_dir,
                               "maize_JIT_edge_all_genes_FCcutoff2_withlog2FC.tsv"))
