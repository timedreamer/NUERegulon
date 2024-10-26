# This script identifies root NxTime genes from maize time-course data.

# Author: Ji Huang
# Date: 2021-04-06

# 0.Prep ----------------------------------------------------------------------------

library(here)
library(edgeR)
library(tidyverse)
library(splines)
library(edge)

source(here("src", "functions", "load_read-count_design-table.R"))

# 1. Prep design matrix and filter genes. -------------------------------------------

Ntreatment <- lib_table$Ntreatment
lib_table$time <- as.integer(as.character(lib_table$time))
X <- ns(lib_table$time, df=5)
design <- model.matrix(~X*Ntreatment)
dim(design)
colnames(design)

dge <- DGEList(counts=raw_data)
keep <- filterByExpr(dge, design)
sum(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# voom transformation.
v <- voom(dge, design, plot=TRUE, normalize="quantile")

vexpr <- v$E
treatment <- lib_table$Ntreatment
time <- lib_table$time


# 2. Call NxTime genes using edge ---------------------------------------------------

de_obj <- build_study(data = vexpr, grp = treatment,
                      tme = time, sampling = "timecourse", basis.df = 5)

de_lrt <- lrt(de_obj, nullDistn = "bootstrap", mod.F = T, seed = 123, bs.its = 1000)
summary(de_lrt)

# Extract gene lists.
sig_results <- qvalueObj(de_lrt)
qvalues <- sig_results$qvalues
gene <- names(sig_results$stat)
cutoff <- 0.05
sigGenes <- qvalues < cutoff
edge_qvalues <- qvalues[sigGenes]
edge_gene <- as.character(gene[sigGenes])
length(edge_gene)

result_df <- tibble(gene = edge_gene, qvalue = edge_qvalues) %>% 
    arrange(qvalue) %>% filter(qvalue < 1E-5)


# 3. Save NxTime result. -------------------------------------------------------

write_tsv(result_df, here("result", "DEgenes", "edge", 
                          "zma_root_NxTgenes_df5_edge.tsv"))

