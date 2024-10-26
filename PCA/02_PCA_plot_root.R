# This script plots PCA plot for root data.

# 0. Prep -----------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(here)
library(PCAtools)
library(scales)
library(dendextend)
library(RColorBrewer)

source(here("src", "functions", "load_read-count_design-table.R"))

# 1. Plot the PCA for all sample. -----------------------------------------

dds_mts <- DESeqDataSetFromMatrix(countData = raw_data,
                                  colData = lib_table,
                                  design = ~ time + Ntreatment)

## normalize sequencing depth.
dds_mts <- estimateSizeFactors(dds_mts)

## remove lowly expressed genes.
keep <- rowSums(counts(dds_mts, normalized = TRUE) != 0) >= 30
sum(keep)
dds_mts <- dds_mts[keep,]


vst_mts<- vst(dds_mts, blind = FALSE)

ntop <-  1500

rv <- rowVars(assay(vst_mts))
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
mat <- assay(vst_mts)[select,]

lib_table <- as.data.frame(lib_table)
rownames(lib_table) <- colnames(mat)

p <- pca(mat, metadata = lib_table)

p1 <- biplot(p, showLoadings = FALSE, lab = NULL, 
             colby = 'time', shape = "Ntreatment",
             legendPosition = 'right',
             title = "MTS Shoot PCA", 
             caption = paste0("Plotted on ", Sys.Date())
             ) +
    scale_colour_viridis_d(direction = -1) +
    coord_fixed()

p4 <- eigencorplot(p,col = brewer_pal(palette = "RdBu")(10), 
                   colCorval = "white", fontCorval = "bold",
                   corFUN = "pearson",
                   corMultipleTestCorrection = "BH",
                   metavars = c("Ntreatment",'time','rep'))

## Save as pdf
pdf(here("result", "pca", paste0("mts_root_pca_all_PCAtools_", 
                                 Sys.Date(),".pdf")
), width = 7, height = 5)
p1

dev.off()

pdf(here("result", "pca", paste0("mts_root_pca_all_PCAtools_eigencorplot_", 
                                 Sys.Date(),".pdf")
), width = 8, height = 3.5)
p4
dev.off()