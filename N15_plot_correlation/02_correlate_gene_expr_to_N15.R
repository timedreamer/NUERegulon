# This script calculates the PCC between shoot and root NxTime gene expression
# with the N15 measurement. It also plots the correlation.

# Author: Ji Huang


# Output:
# shoot_NxTgene_corr_N15.tsv and root_NxTgene_corr_N15.tsv
# shoot_gene_N15_corr_plots_Npath.pdf and root_gene_N15_corr_plots_Npath.pdf

# 0. Prep -----------------------------------------------------------------

library(here)
library(tidyverse)
library(broom)
library(ggforce)

outdir <- here("result", "N15")

source(here("src", "functions", "go_enrich.R"))

## Load all result gene lists.
source(here("src", "functions", "load_all_gene-list.R"))

calc_n15_corr_shoot <- function(i) {
    
    geneid <- rownames(mts_n_shoot)[i]
    cor_result <- cor.test(mts_n_shoot[i,], n15_shoot[1,])
    cor_result <- broom::tidy(cor_result) %>% 
        mutate(geneid = geneid, generow = i)
    
    return(cor_result)
    
}

calc_n15_corr_root <- function(i) {
    
    geneid <- rownames(mts_n_root)[i]
    cor_result <- cor.test(mts_n_root[i,], n15_root[1,])
    cor_result <- broom::tidy(cor_result) %>% 
        mutate(geneid = geneid, generow = i)
    
    return(cor_result)
    
}

# 1. Load N15 -------------------------------------------------------------

n15 <- readxl::read_excel(here("data", "N15", "jh_organized_result.xlsx"),
                  skip = 1,
                  col_names = c("sample_id", "N15", "sample_name", "time", "tissue")) %>% 
    mutate(time_linear = str_remove(time, "m")) %>% 
    mutate(time_linear = as.numeric(time_linear)) %>% 
    mutate(time = factor(time, 
                         levels = c("0m","5m","10m","15m","20m","30m","45m","60m","90m","120m"))) 

## 1.1 N15 shoot
n15_shoot <- n15 %>% 
    filter(tissue == "shoot") %>% 
    dplyr::select(sample_name, N15) %>% 
    mutate(sample_name = str_remove(sample_name, " \\+N")) %>% 
    pivot_wider(names_from = sample_name, values_from = N15) %>% 
    as.data.frame()

rownames(n15_shoot) <- "N15_shoot"

n15_shoot <- as.matrix(n15_shoot)

## 1.2 N15 root
n15_root <- n15 %>% 
    filter(tissue == "root") %>% 
    dplyr::select(sample_name, N15) %>% 
    mutate(sample_name = str_remove(sample_name, " \\+N")) %>% 
    pivot_wider(names_from = sample_name, values_from = N15) %>% 
    as.data.frame()

rownames(n15_root) <- "N15_root"

n15_root <- as.matrix(n15_root)

# 2. Prepare gene expression -------------------------------------------------

## 2.1 Prepare shoot expression
mts_quant <- read_tsv(here("result", "quant_normalized_expr_all.tsv.gz"))

mts_n_shoot <- mts_quant %>% filter(geneid %in% zma_shoot$gene) %>% 
    dplyr::rename(`0mS1-N` = `0mS1`, `0mS2-N` = `0mS2`, `0mS3-N` = `0mS3`) %>%
    dplyr::select(geneid, starts_with("0mS"), ends_with("-N", ignore.case = FALSE)) %>% 
    pivot_longer(col=-geneid, names_to = "condition", values_to = "quant") %>% 
    mutate(condition = paste0("MTS", condition)) %>% 
    tidyr::separate(condition, into = c("condition", "Ntreatment"), sep="-") %>% 
    dplyr::select(-Ntreatment) %>% 
    pivot_wider(names_from = condition, values_from = quant) %>% 
    mutate(MTS120mS3 = as.integer((MTS120mS1 + MTS120mS2)/2))

col_order <- c("geneid", colnames(n15_shoot))

mts_n_shoot <- mts_n_shoot[col_order] %>% 
    as.data.frame()

rownames(mts_n_shoot) <- mts_n_shoot$geneid
mts_n_shoot$geneid <- NULL
mts_n_shoot <- as.matrix(mts_n_shoot)

## 2.2 Prepare root expression

mts_quant <- read_tsv(here("result", 
                           "quant_normalized_expr_all_mts_root.tsv.gz"))

mts_n_root <- mts_quant %>% filter(geneid %in% zma_root$gene) %>% 
    pivot_longer(col=-geneid, names_to = "condition", values_to = "quant") %>% 
    tidyr::separate(condition, into = c("condition", "Ntreatment"), sep="-") %>% 
    tidyr::separate(condition, sep = "mR", into = c("time", "rep")) %>% 
    mutate(time = as.integer(time))

mts_n_root <- mts_quant %>% filter(geneid %in% zma_root$gene) %>% 
    dplyr::rename(`0mR1-N` = `0mR1-noN`, `0mR2-N` = `0mR2-noN`, `0mR3-N` = `0mR3-noN`) %>%
    dplyr::select(geneid, ends_with("-N", ignore.case = FALSE)) %>% 
    pivot_longer(col=-geneid, names_to = "condition", values_to = "quant") %>% 
    mutate(condition = paste0("MTS", condition)) %>% 
    tidyr::separate(condition, into = c("condition", "Ntreatment"), sep="-") %>% 
    dplyr::select(-Ntreatment) %>% 
    pivot_wider(names_from = condition, values_from = quant)

col_order <- c("geneid", colnames(n15_root))

mts_n_root <- mts_n_root[col_order] %>% 
    as.data.frame()

rownames(mts_n_root) <- mts_n_root$geneid
mts_n_root$geneid <- NULL
mts_n_root <- as.matrix(mts_n_root)

# 3. Correlation in Shoot -------------------------------------------------

shoot_corr_all <- map_dfr(1:nrow(mts_n_shoot), 
                          calc_n15_corr_shoot) %>% 
    mutate(padj = p.adjust(.$p.value, method = "BH"))

shoot_gene_anno <- get_geneAnnotation(geneList = shoot_corr_all$geneid, 
                                     species = "zma")

shoot_corr_all <- shoot_corr_all %>% 
    left_join(shoot_gene_anno, by = c("geneid" = "ensembl_gene_id")) %>% 
    dplyr::rename(pcc = estimate)
    
# 4. Correlation in Root --------------------------------------------------

root_corr_all <- map_dfr(1:nrow(mts_n_root), 
                          calc_n15_corr_root) %>% 
    mutate(padj = p.adjust(.$p.value, method = "BH"))

root_gene_anno <- get_geneAnnotation(geneList = root_corr_all$geneid, 
                                     species = "zma")

root_corr_all <- root_corr_all %>% 
    left_join(root_gene_anno, by = c("geneid" = "ensembl_gene_id")) %>% 
    dplyr::rename(pcc = estimate)

root_corr_all %>% filter(padj < 0.01) %>% dim() # 735

# 5. Save output ----------------------------------------------------------

write_tsv(shoot_corr_all, here(outdir, "shoot_NxTgene_corr_N15.tsv"))
write_tsv(root_corr_all, here(outdir, "root_NxTgene_corr_N15.tsv"))


# 6. Plot the 15N-gene expression correlation for Npathway genes -----------------------------------------------------------

shoot_corr_all <- read_tsv(here(outdir, "shoot_NxTgene_corr_N15.tsv"))
root_corr_all <- read_tsv(here(outdir, "root_NxTgene_corr_N15.tsv"))


## Nmetabolism genes from `008_Npathway_gene_combine_with_name.R`.
Npath_genes <- read_tsv(here("result",
                             "zma_Npathway_genes_with_names_20230419.tsv")) %>% 
    distinct()

STime <- paste0("T", c(0,5,10,15,20,30,45,60,90,120))
STime <- rep(STime,each = 3)
STime <- factor(STime, levels = c("T0","T5","T10","T15","T20","T30","T45",
                                  "T60","T90","T120"))

## 8.1 Shoot
shoot_corr_Npath <- Npath_genes |> 
    left_join(shoot_corr_all, by = c("sgene" = "geneid")) |> 
    filter(padj < 0.01) |> arrange(desc(pcc)) |> distinct(sgene, .keep_all = TRUE)


sig_shoot <- shoot_corr_Npath

pdf(here(outdir, "shoot_gene_N15_corr_plots_Npath.pdf"), height = 4, width = 5)

for (i in 1: nrow(sig_shoot)) {
    
    sig_gene <- sig_shoot %>% 
        dplyr::slice(i)
    
    generow = sig_gene$generow
    
    p_title = paste0(sig_gene$sgene, " ", sig_gene$zma_gene_symbol)
    p_subtitle = paste0("PCC: ", round(sig_gene$pcc,2), ". padj: ",
                        format(sig_gene$padj, scientific = TRUE, digits = 2)
                        )
    p_caption = paste0(sig_gene$pathway_step, "\n", sig_gene$description)
    
    
    df1 <- tibble(gene_expr = mts_n_shoot[generow,], shoot_N15 = n15_shoot[1,], 
                  sample = STime)
    
    p1 <- ggplot(df1, aes(gene_expr, shoot_N15)) +
        geom_point(aes(color = STime), size = 2) +
        geom_smooth(se = FALSE, color = "lightblue") +
        scale_colour_viridis_d(direction = -1) +
        xlab("Norm gene expression") + 
        ylab("shoot N15 (%)") +
        labs(title = p_title, subtitle = p_subtitle,
             caption = p_caption) +
        cowplot::theme_cowplot()
    
    plot(p1)
    
}

dev.off()


## 8.2 Root
root_corr_Npath <- Npath_genes |> 
    left_join(root_corr_all, by = c("sgene" = "geneid")) |> 
    filter(padj < 0.01) |> arrange(desc(pcc)) |> distinct(sgene, .keep_all = TRUE)


sig_root <- root_corr_Npath

pdf(here(outdir, "root_gene_N15_corr_plots_Npath.pdf"), height = 4, width = 5)

for (i in 1: nrow(sig_root)) {
    
    sig_gene <- sig_root %>% 
        dplyr::slice(i)
    
    generow = sig_gene$generow
    
    p_title = paste0(sig_gene$sgene, " ", sig_gene$zma_gene_symbol)
    p_subtitle = paste0("PCC: ", round(sig_gene$pcc,2), ". padj: ",
                        format(sig_gene$padj, scientific = TRUE, digits = 2)
                        )
    p_caption = paste0(sig_gene$pathway_step, "\n", sig_gene$description)
    
    
    df1 <- tibble(gene_expr = mts_n_root[generow,], root_N15 = n15_root[1,], 
                  sample = STime)
    
    p1 <- ggplot(df1, aes(gene_expr, root_N15)) +
        geom_point(aes(color = STime), size = 2) +
        geom_smooth(se = FALSE, color = "lightblue")+
        scale_colour_viridis_d(direction = -1) +
        xlab("Norm gene expression") + 
        ylab("root N15 (%)") +
        labs(title = p_title, subtitle = p_subtitle,
             caption = p_caption)+
        cowplot::theme_cowplot()
    
    plot(p1)
    
}

dev.off()

